#!/usr/bin/env python3
"""
Script to discover available metadata values in CELLxGENE Census for breast tissue.
This helps us determine the correct filter values for querying.
"""

import cellxgene_census
import pandas as pd

print("Opening CELLxGENE Census...")
with cellxgene_census.open_soma(census_version="latest") as census:

    # Query for disease values in breast tissue using tissue_general
    print("\n" + "="*60)
    print("Discovering metadata for tissue_general == 'breast'...")
    print("="*60)

    try:
        # Query breast tissue to see available disease values
        disease_df = cellxgene_census.get_obs(
            census,
            "Homo sapiens",
            value_filter="tissue_general == 'breast'",
            column_names=["disease", "tissue", "tissue_general", "assay", "is_primary_data"],
        )

        print("\nUnique 'disease' values in breast tissue:")
        unique_diseases = disease_df['disease'].unique()
        for disease in sorted(unique_diseases):
            count = (disease_df['disease'] == disease).sum()
            print(f"  - '{disease}' ({count:,} cells)")

        print("\nUnique 'tissue' values:")
        unique_tissues = disease_df['tissue'].unique()
        for tissue in sorted(unique_tissues)[:20]:  # Show first 20
            count = (disease_df['tissue'] == tissue).sum()
            print(f"  - '{tissue}' ({count:,} cells)")

        print("\nUnique 'assay' values in breast tissue:")
        unique_assays = disease_df['assay'].unique()
        for assay in sorted(unique_assays):
            count = (disease_df['assay'] == assay).sum()
            print(f"  - '{assay}' ({count:,} cells)")

        # Check specific combinations
        print("\n" + "="*60)
        print("Cell counts for key combinations...")
        print("="*60)

        # Normal + 10x 3' v3 + primary data
        normal_10x = disease_df[(disease_df['disease'] == 'normal') &
                                 (disease_df['assay'] == "10x 3' v3") &
                                 (disease_df['is_primary_data'] == True)]
        print(f"\nNormal breast + 10x 3' v3 + primary: {len(normal_10x):,} cells")

        # Breast cancer variants + 10x 3' v3 + primary data
        cancer_keywords = ['carcinoma', 'cancer', 'adenocarcinoma']
        for keyword in cancer_keywords:
            cancer_10x = disease_df[(disease_df['disease'].str.contains(keyword, case=False, na=False)) &
                                     (disease_df['assay'] == "10x 3' v3") &
                                     (disease_df['is_primary_data'] == True)]
            if len(cancer_10x) > 0:
                diseases_found = cancer_10x['disease'].unique()
                print(f"\n'{keyword}' + 10x 3' v3 + primary: {len(cancer_10x):,} cells")
                print(f"  Disease values: {diseases_found}")

        print("\n" + "="*60)
        print("Discovery complete!")
        print("="*60)

    except Exception as e:
        print(f"\nError: {e}")
        print("\nTrying alternative: tissue_general == 'mammary gland'...")

        try:
            alt_df = cellxgene_census.get_obs(
                census,
                "Homo sapiens",
                value_filter="tissue_general == 'mammary gland'",
                column_names=["disease", "tissue", "assay"],
            )
            print(f"Found {len(alt_df):,} cells with tissue_general == 'mammary gland'")
            print(f"Diseases: {alt_df['disease'].unique()}")
        except Exception as e2:
            print(f"Also failed: {e2}")
