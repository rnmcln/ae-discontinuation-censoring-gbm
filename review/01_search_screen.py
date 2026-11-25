"""Search PubMed Central for single-arm phase II trials in high-grade glioma and
extract candidate passages on AE-related discontinuation and survival follow-up
to support manual coding. Screening decisions are not made here: they were made
manually for every record and are recorded in review_records.csv; the coded data
and verbatim evidence are in review_trials.csv and review_evidence.csv.
NCBI_API_KEY is optional (raises the rate limit). Writes candidates.csv and
screen_auto.csv; full-text XML files are cached in pmc_xml/. Some articles have
no body text in the E-utilities XML; these were read on the PubMed Central web
page (text_source in review_trials.csv).
"""
import json, os, re, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
import pandas as pd

KEY = os.environ.get("NCBI_API_KEY", "")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TERM = ('(glioblastoma[ti] OR "high-grade glioma"[ti] OR gliomas[ti] OR glioma[ti]) AND '
        '("phase II"[ti] OR "phase 2"[ti] OR "phase IIa"[ti] OR "phase Ib/II"[ti] OR "phase I/II"[ti]) AND '
        '("single-arm" OR "single arm" OR "open-label, non-randomised" OR "nonrandomized" OR "non-randomized") '
        'AND "overall survival" AND 2015:2025[pdat]')
TITLE_FLAG = re.compile(  # flags titles for manual screening; not an exclusion rule
   
    r"pediatric|paediatric|child|DIPG|pontine|protocol|randomi[sz]ed|meta-analysis|review|low-grade|low grade|"
    r"retrospective|secondary analysis|post hoc|exploratory analysis|imaging|MRI|biomarker|radiomic|"
    r"Grade II Glioma|Regulatory Science|Inability of positive|ACRIN|Iron-based|Long-Term Report", re.I)
PATTERNS = {
    "ae_disc": r"(discontinu|stopp|withdr|terminat)\w*[^.]{0,120}(adverse|toxicit|side effect|AE\b|AEs\b|intoleran)"
               r"|(adverse|toxicit|AE\b|AEs\b|intoleran)\w*[^.]{0,120}(discontinu|withdr|stopp|terminat)",
    "ltfu": r"lost to follow|loss to follow|withdr\w* (of )?consent|consent withdraw|withdrew consent",
    "surv_fu": r"(follow\w*[^.]{0,60}(for|until) (survival|death))|(survival (status|follow))"
               r"|(after (treatment |study )?discontinu\w*[^.]{0,80}follow)|(follow\w*[^.]{0,80}after (treatment )?discontinu)",
    "censor": r"censor",
    "ic_sens": r"informative censoring|sensitivity analys|inverse probability|competing risk",
    "primary_endpoint": r"primary (?:end ?point|objective|outcome|efficacy end ?point|aim)",
    "comparator": r"historical|benchmark|external control|matched control|synthetic control|control cohort",
}


def get(url):
    for _ in range(3):
        try:
            return urllib.request.urlopen(url, timeout=60).read()
        except Exception:
            time.sleep(2)
    raise RuntimeError(url)


def esearch(term, n=400):
    q = dict(db="pmc", term=term, retmax=n, retmode="json")
    if KEY:
        q["api_key"] = KEY
    return json.loads(get(EU + "esearch.fcgi?" + urllib.parse.urlencode(q)))["esearchresult"]["idlist"]


def efetch(pid):
    os.makedirs("pmc_xml", exist_ok=True)
    fn = f"pmc_xml/{pid}.xml"
    if os.path.exists(fn):
        return open(fn).read()
    q = dict(db="pmc", id=pid, retmode="xml")
    if KEY:
        q["api_key"] = KEY
    t = get(EU + "efetch.fcgi?" + urllib.parse.urlencode(q)).decode()
    open(fn, "w").write(t)
    time.sleep(0.12)
    return t


def _clean(s):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", s))


def body_text(xml):
    """Abstract and body text, tables removed."""
    a = re.search(r"<abstract.*?</abstract>", xml, re.S)
    m = re.search(r"<body.*?</body>", xml, re.S)
    s = (a.group(0) if a else "") + " " + (m.group(0) if m else xml)
    s = re.sub(r"<table-wrap.*?</table-wrap>", " ", s, flags=re.S)
    return _clean(s)


def table_text(xml):
    """Tables, table footnotes and figure captions."""
    parts = re.findall(r"<table-wrap.*?</table-wrap>", xml, re.S) + re.findall(r"<caption.*?</caption>", xml, re.S)
    return _clean(" ".join(parts))


def main():
    rows = []
    for pid in esearch(TERM):
        xml = efetch(pid)
        try:
            r = ET.fromstring(xml)
        except ET.ParseError:
            continue
        art = r.find(".//article")
        title = "".join(r.find(".//article-title").itertext()) if r.find(".//article-title") is not None else ""
        ids = {x.get("pub-id-type"): x.text for x in r.findall(".//article-id")}
        rows.append(dict(pmcid="PMC" + pid, pmid=ids.get("pmid", ""), doi=ids.get("doi", ""),
                         year=r.findtext(".//pub-date/year") or "", journal=r.findtext(".//journal-title") or "",
                         type=art.get("article-type") if art is not None else "",
                         title=re.sub(r"\s+", " ", title)))
    cand = pd.DataFrame(rows)
    NON_REPORT = {"abstract", "meeting-report", "editorial", "letter", "reply", "review-article", "systematic-review"}
    cand["title_flag"] = cand.title.fillna("").str.contains(TITLE_FLAG)
    cand.to_csv("candidates.csv", index=False)
    keep = cand[~cand.type.isin(NON_REPORT)].copy()
    out = []
    for _, r in keep.iterrows():
        xml = efetch(r.pmcid[3:])
        sents = re.split(r"(?<=[.;])\s+", body_text(xml))
        tsents = re.split(r"(?<=[.;|])\s+", table_text(xml))
        rec = r.to_dict()
        for k, p in PATTERNS.items():
            hits = [s for s in sents if re.search(p, s, re.I)]
            thits = [s for s in tsents if re.search(p, s, re.I)]
            rec[k + "_n"] = len(hits)
            rec[k + "_snippets"] = " || ".join(h[:300] for h in hits[:4])
            rec[k + "_table_n"] = len(thits)
            rec[k + "_table_snippets"] = " || ".join(h[:300] for h in thits[:3])
        out.append(rec)
    pd.DataFrame(out).to_csv("screen_auto.csv", index=False)
    print(f"{len(cand)} records retrieved; candidate passages extracted for {len(out)} full reports")


if __name__ == "__main__":
    main()
