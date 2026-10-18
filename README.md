# MigraTarget

Exploratory tool that generates **drug-target hypotheses for migraine** from pharmacovigilance data.
Built as a learning project in Python + Streamlit, with every design decision and limitation documented.

> **Not medical advice.** Outputs are statistical signals from spontaneous-report data. They do not show that a drug causes or treats migraine.

## What it does

**Phase 1 - FAERS signal detection**

1. Downloads migraine adverse-event reports from the FDA's FAERS database through the openFDA API (sample: the 500 most recent reports).
2. Keeps reports where the drug is the *primary suspect* and removes duplicate reports.
3. Counts drugs and keeps those that appear in at least 10 reports.
4. Removes non-pharmacological entries (vaccines, sanitisers).
5. Computes the reporting odds ratio (ROR) with a 95% confidence interval for each candidate, against the **full** FAERS database.
6. Keeps statistically significant signals and excludes drugs with already-known migraine mechanisms (CGRP, triptans, analgesics, etc.).
7. Classifies the remaining candidates by ATC class (RxNorm/RxClass) and charts the distribution.
8. Saves an audit table with the reason each candidate was kept or discarded.

**Phase 2 - Literature check**
For each final candidate, searches PubMed (`(drug) AND (Migraine*)`) and lists the most recent articles with date, authors and link.

**Chat**
Ask in natural language (e.g. "What is the ROR of methotrexate?"). Gemini decides which tool to call (ROR calculation or PubMed search) and answers with live data.


## Limitations

- Candidates are discovered from a 500-report sample, not the full database (the ROR itself is always computed against the full database).
- The list of "known-mechanism" drugs is curated by hand and not exhaustive.
- ATC classification uses the first class returned by the API, which is not always the most representative one.
- The PubMed check reads only the 10 most recent articles per drug, and brand names may be searched less effectively than generic names.
- The chat uses a free Gemini tier with a daily quota.
- Disproportionality analysis cannot establish causality.

## Setup

1. Install Python 3 and the dependencies: `pip install -r requirements.txt`
2. Copy `.streamlit/secrets.toml.example` to `.streamlit/secrets.toml` and add your own keys:
   - openFDA: https://open.fda.gov/apis/authentication/
   - NCBI (PubMed): account settings > API Key Management
   - Google AI Studio: https://aistudio.google.com/apikey
3. Run: `streamlit run app.py`

## Built with

Python, pandas, Streamlit, openFDA, PubMed (Biopython), RxNorm/RxClass, Gemini.
