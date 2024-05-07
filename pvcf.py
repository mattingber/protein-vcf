import json
import click
import re
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback

LETTER_REGEX = r"^(\d+)(.*)"
INSERT_OPERATOR = ">"
DEL_OPERATOR = "del"
DEL_STRING = "<DEL>"


class Diff(BaseModel):
    pos: str
    ref: str
    alt: Optional[str]
    info: Optional[str]


class Haplotype(BaseModel):
    frequency: float
    samples: Dict[str, int]
    aligned_sequences: List[str]
    diffs: List[Diff]
    name: str


class Item(BaseModel):
    prot: str
    pos: str
    ref: str
    alt: List[str]
    info: str
    samples: Dict[str, List[int]]


HEADER_TO_KEY_MAP = {
    'PROT': 'prot',
    'POS': 'pos',
    'REF': 'ref',
    'ALT': 'alt',
    'INFO': 'info'
}

ALT_SEPARATOR = ','
SAMPLE_SEPARATOR = '|'
TAB_STRING = "\t"

def get_formatted_haplotypes(json_content: dict) -> tuple[List[Haplotype], List[str]]:
    sample_ids = set()
    protein_haplotypes = json_content.get("protein_haplotypes", [])
    formatted_haplotypes = []
    for haplotype in protein_haplotypes:
        frequency = haplotype.get("frequency", 0)
        samples = haplotype.get("samples", {})
        sample_ids.update(samples.keys())

        diffs = []
        aligned_sequences = haplotype.get("aligned_sequences", [""])
        for diff_data in haplotype.get("diffs", []):
            diff = diff_data.get("diff")
            formatted_diff = format_diff(diff, aligned_sequences[0])
            if not formatted_diff:
                continue
            diffs.append(formatted_diff)

        name = haplotype.get("name", "")
        formatted_haplotype = Haplotype(
            frequency=frequency,
            samples=samples,
            aligned_sequences=aligned_sequences,
            diffs=diffs,
            name=name,
        )
        formatted_haplotypes.append(formatted_haplotype)

    return formatted_haplotypes, list(sample_ids)


def format_diff(diff: str | None, ref_sequence: str) -> Diff | None:
    if not diff:
        return None

    matches = re.match(LETTER_REGEX, diff)
    if matches:
        pos = matches.group(1)
        remaining_part = matches.group(2)
        if INSERT_OPERATOR in remaining_part:
            [ref, alt] = remaining_part.split(INSERT_OPERATOR)
            return Diff(pos=pos, ref=ref, alt=alt, info="")

        if DEL_OPERATOR in remaining_part:
            [_, deletion_part] = remaining_part.split(DEL_OPERATOR)
            ref = ref_sequence[int(pos)]
            alt = DEL_STRING
            deletion_part_range = deletion_part.strip('{').strip('}')
            try:
                deletion_part_range_int = int(deletion_part_range)
                info = (
                    f"SVTYPE=DEL;END={int(pos) + deletion_part_range_int}"
                )
            except ValueError:
                info = (
                    f"SVTYPE=DEL; SVID={pos}del{deletion_part_range}"
                )
            return Diff(pos=pos, ref=ref, alt=alt, info=info)

        return None
    return None


def handle_samples(prev_samples: dict, haplotype: Haplotype, alt_index: int):
    samples: Dict[str, List] = {**prev_samples}
    for sample_id, value in haplotype.samples.items():
        prev_sample = prev_samples.get(sample_id)
        samples_value = prev_sample if prev_sample else []
        for _ in range(value):
            samples_value.append(int(alt_index))
            samples[sample_id] = samples_value
    return samples


def build_items(formatted_haplotypes: List[Haplotype]) -> List[Item]:
    vcf_rows = {}
    for haplotype in formatted_haplotypes:
        for diff in haplotype.diffs:
            row_key = f"{diff.pos}:{diff.ref}"
            prev_row = dict(vcf_rows.get(row_key, {}))

            prot = haplotype.name.split(":")[0]
            pos = diff.pos
            ref = diff.ref
            prev_info: str = prev_row.get("info", "")
            if prev_info:
                info = prev_info
            else:
                info = diff.info
            current_alt = diff.alt

            alt: List[str] = prev_row.get("alt", [])
            if current_alt not in alt:
                alt.append(current_alt)

            prev_samples = prev_row.get("samples", {})
            samples: Dict[str, List] = handle_samples(
                prev_samples, haplotype, alt_index=len(alt)
            )
            vcf_rows[row_key] = Item(
                pos=pos,
                prot=prot,
                ref=ref,
                alt=alt,
                info=info,
                samples=samples,
            )
    return vcf_rows.values()


def append_samples_to_items(
    rows: List[Item], all_sample_ids: List[str]
) -> List[Item]:
    rows_with_samples = []
    for row in rows:
        row.samples = {
            **{sample_id: [0, 0] for sample_id in all_sample_ids},
            **row.samples,
        }
        for sample_id, value in row.samples.items():
            if len(value) < 2:
                row.samples[sample_id].extend([0] * (2 - len(value)))
        rows_with_samples.append(row)
    return rows_with_samples

def generate_vcf_rows(items: List[Item], headers: List[str]):
    rows = []
    for item in items:
        row = []
        dict_item = dict(item)
        for header in headers:
            key = HEADER_TO_KEY_MAP.get(header)
            if key:
                cell = dict_item[key]
                if header == 'ALT':
                    cell = ALT_SEPARATOR.join(cell)
            else:
                sample = item.samples.get(header)
                cell = SAMPLE_SEPARATOR.join(str(num) for num in sample)

            row.append(cell)
        rows.append(row)
    return rows


def write_vcf_row(file, row):
    row_string = TAB_STRING.join(row)
    file.write(row_string + "\n")

def build_vcf_file(vcf_table: List[List[str]], output: str):
    with open(output, "a") as file:
        for row in vcf_table:
            write_vcf_row(file, row)



@click.group("pvcf")
def pvcf():
    """An all-in-one tool for all pvcf related operations."""
    pass


@pvcf.command("convert")
@click.option("--path", "path", "-p", required=True)
@click.option("--output", "output", "-o", required=True)
def convert_json_to_pvcf(path: str, output: str):
    """Convert Haplosaurus JSON output to VCF format."""
    vcf_headers = []
    with open(path, "r") as file:
        try:
            print(f"Converting {path} to VCF format")
            all_transcripts_data = json.load(file)
            for transcript_data in all_transcripts_data:
                formatted_haplotypes, sample_ids = get_formatted_haplotypes(transcript_data)
                items = build_items(formatted_haplotypes)
                items_with_samples = append_samples_to_items(items, sample_ids)
                if not vcf_headers:
                    vcf_headers = [*HEADER_TO_KEY_MAP.keys(), *[str(sample_id) for sample_id in sample_ids]]
                    with open(output, "a") as output_file:
                        write_vcf_row(output_file, vcf_headers)
                vcf_rows = generate_vcf_rows(items_with_samples, vcf_headers)
                # vcf_table = [vcf_headers, *vcf_rows]
                build_vcf_file(vcf_rows, output)

        except Exception as e:
            print(e)
            print(traceback.format_exc())
            raise click.ClickException("Error converting JSON to VCF")


if __name__ == "__main__":
    pvcf()
