"""Abstract-level PubMed query used to judge how representative the open-access
full-text sample is. Writes pubmed_denominator.csv."""
import json, os, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
import pandas as pd

KEY = os.environ.get("NCBI_API_KEY", "")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TERM = ('(glioblastoma[ti] OR "high-grade glioma"[ti] OR glioma[ti] OR gliomas[ti]) AND '
        '("phase II"[ti] OR "phase 2"[ti] OR "phase IIa"[ti] OR "phase Ib/II"[ti] OR "phase I/II"[ti]) AND '
        '2015:2025[dp] AND ("single-arm"[tiab] OR "single arm"[tiab] OR nonrandomized[tiab] OR '
        '"non-randomized"[tiab] OR "non-randomised"[tiab] OR "open-label"[tiab]) NOT '
        '(randomized controlled trial[pt] OR protocol[ti] OR pediatric[ti] OR children[ti])')


def q(path, **kw):
    if KEY:
        kw["api_key"] = KEY
    return urllib.request.urlopen(EU + path + "?" + urllib.parse.urlencode(kw), timeout=60).read()


ids = json.loads(q("esearch.fcgi", db="pubmed", term=TERM, retmax=500, retmode="json"))["esearchresult"]["idlist"]
root = ET.fromstring(q("efetch.fcgi", db="pubmed", id=",".join(ids), retmode="xml"))
rows = []
for a in root.findall(".//PubmedArticle"):
    aid = {x.get("IdType"): x.text for x in a.findall(".//ArticleIdList/ArticleId")}
    rows.append(dict(pmid=a.findtext(".//PMID"), pmcid=aid.get("pmc", ""), journal=a.findtext(".//ISOAbbreviation"),
                     year=a.findtext(".//JournalIssue/PubDate/Year"), title="".join(a.find(".//ArticleTitle").itertext())))
d = pd.DataFrame(rows)
d.to_csv("pubmed_denominator.csv", index=False)
print(f"{len(d)} records; {(d.pmcid != '').sum()} with PubMed Central full text")
print(d[d.pmcid == ""].journal.value_counts().to_string())
