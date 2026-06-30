#!/usr/bin/env python3

import csv
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def fetch_biosample_xml(biosample):
    """Search NCBI BioSample and return the record as XML text."""

    search = subprocess.run(
        ["esearch", "-db", "biosample", "-query", biosample],
        capture_output=True,
        text=True,
        check=True,
    )

    fetch = subprocess.run(
        ["efetch", "-format", "xml"],
        input=search.stdout,
        capture_output=True,
        text=True,
        check=True,
    )

    return fetch.stdout


def main():
    if len(sys.argv) != 2:
        print("Usage: python add_biosample_attributes.py runinfo.csv")
        sys.exit(1)

    runinfo_csv = Path(sys.argv[1])
    output_tsv = runinfo_csv.with_suffix("").with_name(
        runinfo_csv.stem + ".biosample_attributes.tsv"
    )

    # A set prevents duplicate Run/BioSample combinations.
    run_biosamples = set()

    with open(runinfo_csv, newline="") as infile:
        reader = csv.DictReader(infile)

        for row in reader:
            run = row.get("Run", "").strip()
            biosample = row.get("BioSample", "").strip()

            if run and biosample:
                run_biosamples.add((run, biosample))

    with open(output_tsv, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Run", "BioSample", "Attribute", "Value"])

        for run, biosample in sorted(run_biosamples):
            print(f"Fetching BioSample attributes for {run} / {biosample}")

            try:
                xml_text = fetch_biosample_xml(biosample)
                root = ET.fromstring(xml_text)

                for element in root.iter():
                    # Handles XML tags whether or not they include a namespace.
                    tag = element.tag.split("}")[-1]

                    if tag != "Attribute":
                        continue

                    attribute = element.attrib.get("attribute_name", "")
                    value = (element.text or "").strip()

                    if attribute or value:
                        writer.writerow([run, biosample, attribute, value])

            except subprocess.CalledProcessError:
                print(f"Warning: could not fetch {biosample}", file=sys.stderr)
            except ET.ParseError:
                print(f"Warning: invalid XML returned for {biosample}", file=sys.stderr)

    print(f"Wrote: {output_tsv}")


if __name__ == "__main__":
    main()