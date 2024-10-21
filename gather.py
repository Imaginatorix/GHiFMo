from chembl_webresource_client.new_client import new_client
import csv
import os

def last_entry(location):
    # Log location
    location = location+".txt"

    # Create if it does not exist
    if not os.path.exists(location):
        with open(location, "w+") as f:
            f.write("0")
        return 0

    with open(location, "r") as f:
        return int(f.readline())


def _or(a, b, *args):
    try:
        for arg in args:
            a = a[arg]
        return a
    except Exception:
        return b


# Save to CSV
def write2csv(location, batch):
    with open(location, "a+", newline='') as f:
        csv_writer = csv.writer(f)

        start = last_entry(location) == 0
        for data in batch:
            if start:
                # Add headers
                header = ["chembl_id",
                          "smiles",
                          "alogp",
                          "aromatic_rings",
                          "cx_logd",
                          "cx_logp",
                          "cx_most_apka",
                          "cx_most_bpka",
                          "full_molformula",
                          "full_mwt",
                          "hba",
                          "hba_lipinski",
                          "hbd",
                          "hbd_lipinski",
                          "heavy_atoms",
                          "molecular_species",
                          "mw_freebase",
                          "mw_monoisotopic",
                          "np_likeness_score",
                          "num_lipinski_ro5_violations",
                          "num_ro5_violations",
                          "psa",
                          "qed_weighted",
                          "ro3_pass",
                          "rtb"]
                csv_writer.writerow(header)
                start = False

            row = [ _or(data, None, "molecule_chembl_id"),
                    _or(data, None, "molecule_structures", "canonical_smiles"),
                    _or(data, None, "molecule_properties", "alogp"),
                    _or(data, None, "molecule_properties", "aromatic_rings"),
                    _or(data, None, "molecule_properties", "cx_logd"),
                    _or(data, None, "molecule_properties", "cx_logp"),
                    _or(data, None, "molecule_properties", "cx_most_apka"),
                    _or(data, None, "molecule_properties", "cx_most_bpka"),
                    _or(data, None, "molecule_properties", "full_molformula"),
                    _or(data, None, "molecule_properties", "full_mwt"),
                    _or(data, None, "molecule_properties", "hba"),
                    _or(data, None, "molecule_properties", "hba_lipinski"),
                    _or(data, None, "molecule_properties", "hbd"),
                    _or(data, None, "molecule_properties", "hbd_lipinski"),
                    _or(data, None, "molecule_properties", "heavy_atoms"),
                    _or(data, None, "molecule_properties", "molecular_species"),
                    _or(data, None, "molecule_properties", "mw_freebase"),
                    _or(data, None, "molecule_properties", "mw_monoisotopic"),
                    _or(data, None, "molecule_properties", "np_likeness_score"),
                    _or(data, None, "molecule_properties", "num_lipinski_ro5_violations"),
                    _or(data, None, "molecule_properties", "num_ro5_violations"),
                    _or(data, None, "molecule_properties", "psa"),
                    _or(data, None, "molecule_properties", "qed_weighted"),
                    _or(data, None, "molecule_properties", "ro3_pass"),
                    _or(data, None, "molecule_properties", "rtb")]
            csv_writer.writerow(row)


if __name__ == "__main__":
    if not os.path.exists("data"):
        os.mkdir("data")

    # Constants
    DATABASE_LOCATION = "data/molecules.csv"
    BATCH_SIZE = 1000

    # All human molecule-target assays with single target protein assigned
    assays = new_client.assay.filter(assay_organism="Homo sapiens", confidence_score=9).only("assay_chembl_id")

    # Get molecules associated with the assay through the activities
    # All human molecule-target activities of the particular assay
    activities = new_client.activity.filter(target_organism="Homo sapiens").only(["assay_chembl_id", "molecule_chembl_id"])
    # All molecules
    molecules = new_client.molecule.filter(molecule_type="Small molecule").only(["molecule_chembl_id", "molecule_properties", "molecule_structures"])

    # The loop to get the molecules
    last_index = last_entry(DATABASE_LOCATION)
    total_counts = len(assays)

    while last_index < total_counts:
        assay_ids = [assay["assay_chembl_id"] for assay in assays[last_index:last_index+BATCH_SIZE]]
        # Get all molecules associated with the assay
        molecule_ids = list(set([activity["molecule_chembl_id"] for activity in activities.filter(assay_chembl_id__in=assay_ids)]))
        batch = [molecule for molecule in molecules.filter(molecule_chembl_id__in=molecule_ids)]

        # Update last index
        last_index = last_index+BATCH_SIZE
        # Save molecule batch
        write2csv(DATABASE_LOCATION, batch)
        # Log assay number
        with open(DATABASE_LOCATION+".txt", "w") as f:
            f.write(str(last_index))
        print(f"{last_index}/{total_counts} assays scanned...")

