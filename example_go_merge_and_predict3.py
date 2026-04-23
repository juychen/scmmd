#!/usr/bin/env python3
"""
GO term mapping and prediction integration script
First perform mapping, then prediction, finally generate new column: use mapping result if available, otherwise use prediction result
"""

import pandas as pd
import os
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import re
import argparse
import joblib
import json


def preprocess_go_term(go_term):
    """Clean GO term by removing prefixes and special characters"""
    # Remove GO:BP_, GO:MF_, GO:CC_ prefixes if present
    clean_term = re.sub(r'^GO:[A-Z]+_', '', go_term)
    # Remove any remaining special characters and make lowercase
    clean_term = re.sub(r'[^\w\s]', '', clean_term).lower()
    return clean_term


def create_label_mapping(labels):
    """
    Create mapping from string labels to numeric labels and vice versa

    Args:
        labels: Series or array of labels (can be string or numeric)

    Returns:
        tuple: (numeric_labels, str_to_num_mapping, num_to_str_mapping)
    """
    unique_labels = pd.Series(labels).unique()

    # Check if labels are already numeric
    if pd.api.types.is_numeric_dtype(labels):
        # Already numeric, create string mapping for consistency
        str_to_num = {str(label): label for label in unique_labels}
        num_to_str = {label: str(label) for label in unique_labels}
        return labels, str_to_num, num_to_str

    # Create string to numeric mapping
    str_to_num = {label: idx for idx, label in enumerate(sorted(unique_labels))}
    num_to_str = {idx: label for label, idx in str_to_num.items()}

    # Convert labels to numeric
    numeric_labels = pd.Series(labels).map(str_to_num)

    print(f"Created label mapping: {len(str_to_num)} unique labels")
    print(f"String to Numeric mapping: {str_to_num}")

    return numeric_labels, str_to_num, num_to_str


def train_single_classifier(df, group_name="All"):
    """
    Train a single classifier for the given dataframe

    Args:
        df: DataFrame with preprocessed data (must have 'clean_go_term' and 'numeric_labels')
        group_name: Name of the group for logging

    Returns:
        Trained pipeline
    """
    # Split into features (X) and targets (y)
    # copy the df 3 times
    df = pd.concat([df.copy() for _ in range(3)], ignore_index=True)
    X = df['clean_go_term']
    y_module = df['numeric_labels']  # Use numeric labels
    print(df['numeric_labels'].value_counts())
    # Split data into training and test sets
    if len(df) < 4:  # Too few samples for proper train/test split
        print(f"Warning: Too few samples ({len(df)}) for {group_name} group. Using all data for training.")
        X_train, X_test = X, X
        y_module_train, y_module_test = y_module, y_module
    else:
        X_train, X_test, y_module_train, y_module_test = train_test_split(
            X, y_module, test_size=0.01, random_state=42, stratify=y_module
        )

    # Create pipeline for module classification
    module_pipeline = Pipeline([
        ('tfidf', TfidfVectorizer(ngram_range=(1, 2), min_df=1, max_df=0.8, max_features=5000)),
        ('clf', SVC(kernel='linear', C=5, class_weight=None))  # Using SVM for classification
    ])

    # Train classifier
    module_pipeline.fit(X_train, y_module_train)

    # Evaluate performance if we have test data
    if len(X_test) > 0 and not X_train.equals(X_test):
        module_pred = module_pipeline.predict(X_test)
        print(f"{group_name} Module Classification Report:")
        print(classification_report(y_module_test, module_pred))

    return module_pipeline


def train_go_classifier(csv_file_path, input_col=None, output_col=None, neurotransmitter_col="Neurotransmitter"):
    """
    Train separate classifiers for NN and N groups based on neurotransmitter column

    Args:
        csv_file_path: Path to CSV file with columns: go_term, module, function, neurotransmitter
        input_col: Name of the input feature column (GO terms)
        output_col: Name of the output label column (modules/categories)
        neurotransmitter_col: Name of the neurotransmitter column

    Returns:
        Dictionary containing trained pipelines for NN and N groups and label mappings
    """
    # Load and preprocess data
    if csv_file_path.endswith('.xlsx') or csv_file_path.endswith('.xls'):
        df = pd.read_excel(csv_file_path, sheet_name=0)
    else:
        df = pd.read_csv(csv_file_path)
    df = df.dropna()

    # Determine input and output columns
    if input_col is None:
        input_col = df.columns[0]
        print(f"Using first column '{input_col}' as input features")

    if output_col is None:
        output_col = df.columns[1]
        print(f"Using second column '{output_col}' as output labels")

    # Check if columns exist
    if input_col not in df.columns:
        raise ValueError(f"Input column '{input_col}' not found in data")
    if output_col not in df.columns:
        raise ValueError(f"Output column '{output_col}' not found in data")

    # Preprocess input features
    df['clean_go_term'] = df[input_col].apply(preprocess_go_term)

    # Check if neurotransmitter column exists
    if neurotransmitter_col not in df.columns:
        print(f"Warning: {neurotransmitter_col} column not found. Using all data for single model.")
        # Create global label mapping for single model
        numeric_labels, str_to_num_mapping, num_to_str_mapping = create_label_mapping(df[output_col])
        df['numeric_labels'] = numeric_labels
        result = train_single_classifier(df, group_name="All")
        return {
            'classifier': result,
            'str_to_num_mapping': str_to_num_mapping,
            'num_to_str_mapping': num_to_str_mapping
        }

    # Split data based on neurotransmitter type
    df_nn = df[df[neurotransmitter_col] == 'NN'].copy()
    df_n = df[df[neurotransmitter_col] != 'NN'].copy()

    print(f"NN group samples: {len(df_nn)}")
    print(f"N group samples: {len(df_n)}")

    classifiers = {}

    # Train NN model if data exists
    if len(df_nn) > 0:
        print("\n=== Training NN Model ===")
        # Create separate label mapping for NN group
        nn_numeric_labels, nn_str_to_num, nn_num_to_str = create_label_mapping(df_nn[output_col])
        df_nn['numeric_labels'] = nn_numeric_labels
        classifiers['NN'] = train_single_classifier(df_nn, group_name="NN")
        classifiers['NN_str_to_num_mapping'] = nn_str_to_num
        classifiers['NN_num_to_str_mapping'] = nn_num_to_str

    # Train N model if data exists
    if len(df_n) > 0:
        print("\n=== Training N Model ===")
        # Create separate label mapping for N group
        n_numeric_labels, n_str_to_num, n_num_to_str = create_label_mapping(df_n[output_col])
        df_n['numeric_labels'] = n_numeric_labels
        classifiers['N'] = train_single_classifier(df_n, group_name="N")
        classifiers['N_str_to_num_mapping'] = n_str_to_num
        classifiers['N_num_to_str_mapping'] = n_num_to_str

    return classifiers


def merge_go_data(go_dir_path, go_group_path, df_class_path, output_path,
                  duplicate_handling='filter', nlogp_threshold=0.7):
    """
    Merge GO data and add grouping information

    Args:
        go_dir_path: GO enrichment results directory
        go_group_path: GO grouping mapping file path
        df_class_path: Class mapping file path
        output_path: Output file path
        duplicate_handling: How to handle duplicates ('filter', 'nlogp_compare', 'nlogp_keep_highest')
        nlogp_threshold: Threshold for nlogp difference comparison (0.0-1.0, e.g., 0.7 for 70%)

    Returns:
        Merged DataFrame
    """
    print("=== Step 1: Merging GO data ===")

    # Read grouping mapping file
    go_group_df = pd.read_excel(go_group_path, sheet_name=0)

    # Read class mapping file
    df_Class = pd.read_excel(df_class_path)
    subclass_to_class_map = df_Class.drop_duplicates(subset=['Subclass']).set_index('Subclass')[
        'Class_for_degs'].to_dict()

    # Merge all GO files
    go_all = []
    for file in os.listdir(go_dir_path):
        file_path = os.path.join(go_dir_path, file)
        if ('enrichment' in file) and (file.endswith('Down.csv') or file.endswith('Up.csv')):
            df = pd.read_csv(file_path)
            if ('source' not in df.columns) and ('ONTOLOGY' in df.columns ):
                df['source'] = 'GO:'+ df['ONTOLOGY']
            if ('term_name' not in df.columns) and ('Description' in df.columns):
                df['term_name'] = df['Description']
            if ('p_value') not in df.columns and ('pvalue' in df.columns):
                df['p_value'] = df['pvalue']
            go_all.append(df)
    go_all_df = pd.concat(go_all, axis=0, ignore_index=True)
    go_all_df['GO'] = go_all_df['source'] + '_' + go_all_df['term_name']

    # Process N and NN groups separately

    
    go_all_N = go_all_df[go_all_df['Neurotransmitter'] != 'NN'].copy()
    go_all_NN = go_all_df[go_all_df['Neurotransmitter'] == 'NN'].copy()

    if 'Neurotransmitter' not in go_group_df.columns:
        go_group_df_N = go_group_df.copy()
        go_group_df_NN = go_group_df.copy()
    else:
        go_group_df_N = go_group_df[go_group_df['Neurotransmitter'] != 'NN'].copy()
        go_group_df_NN = go_group_df[go_group_df['Neurotransmitter'] == 'NN'].copy()

    # Perform mapping
    go_all_N = go_all_N.merge(go_group_df_N[['GO', 'Group']], on='GO', how='left')
    go_all_NN = go_all_NN.merge(go_group_df_NN[['GO', 'Group']], on='GO', how='left')
    go_all_merged = pd.concat([go_all_N, go_all_NN], axis=0, ignore_index=True)

    # Add Class information
    go_all_merged['Class'] = go_all_merged['Subclass'].map(subclass_to_class_map)
    if 'Region Subclass' in go_all_merged.columns:
        go_all_merged.rename(columns={'Region Subclass': 'Region subclass'}, inplace=True)
    if 'Region' not in go_all_merged.columns:
        go_all_merged['Region'] = go_all_merged['Region subclass'].str.split(' ').str[0]

    # Handle duplicate rows per strategy
    duplicates = go_all_merged.duplicated(subset=["Region subclass", "Sex", "GO"], keep=False)
    duplicate_count = int(duplicates.sum())
    if duplicate_count > 0:
        print(f"Found {duplicate_count} duplicate records")
        if duplicate_handling == 'filter':
            print("Removing all duplicate records (filter strategy)...")
            go_all_merged = go_all_merged[~duplicates].copy()
        elif duplicate_handling in ('Up', 'Down'):
            # Ensure nlogp exists for sign filtering
            if 'nlogp' not in go_all_merged.columns and 'p_value' in go_all_merged.columns:
                go_all_merged['nlogp'] = -np.log10(go_all_merged['p_value'])
                if 'Direction' in go_all_merged.columns:
                    go_all_merged['nlogp'] = go_all_merged.apply(
                        lambda row: -row['nlogp'] if row['Direction'] == 'Down' else row['nlogp'], axis=1
                    )
            if 'nlogp' not in go_all_merged.columns:
                print("Warning: No nlogp/p_value available. Falling back to filter strategy.")
                go_all_merged = go_all_merged[~duplicates].copy()
            else:
                key_cols = ["Region subclass", "Sex", "GO"]
                to_drop = []
                want_positive = (duplicate_handling == 'Up')
                for keys, group in go_all_merged[duplicates].groupby(key_cols):
                    mask = group['nlogp'] > 0 if want_positive else group['nlogp'] < 0
                    chosen = group[mask]
                    if len(chosen) == 0:
                        # If no records match the requested sign, drop all duplicates in this group
                        to_drop.extend(list(group.index))
                    else:
                        # Keep all matching sign; drop the rest in the duplicate group
                        to_drop.extend([idx for idx in group.index if idx not in chosen.index])
                if to_drop:
                    before = len(go_all_merged)
                    go_all_merged = go_all_merged.drop(index=list(set(to_drop)))
                    after = len(go_all_merged)
                    print(f"Duplicates resolved by sign ('{duplicate_handling}'): removed {before - after} rows")
                else:
                    print("Duplicate sign filtering removed 0 rows")
        else:
            # Ensure nlogp exists when needed
            if duplicate_handling in ('nlogp_compare', 'nlogp_keep_highest'):
                if 'nlogp' not in go_all_merged.columns:
                    if 'p_value' in go_all_merged.columns:
                        go_all_merged['nlogp'] = -np.log10(go_all_merged['p_value'])
                        if 'Direction' in go_all_merged.columns:
                            go_all_merged['nlogp'] = go_all_merged.apply(
                                lambda row: -row['nlogp'] if row['Direction'] == 'Down' else row['nlogp'], axis=1
                            )
                    else:
                        print("Warning: No p_value column found. Falling back to filter strategy.")
                        go_all_merged = go_all_merged[~duplicates].copy()
                        duplicates = go_all_merged.duplicated(subset=["Region subclass", "Sex", "GO"], keep=False)

            # Group duplicates and resolve
            key_cols = ["Region subclass", "Sex", "GO"]
            dup_df = go_all_merged.copy()
            to_drop = []
            for keys, group in dup_df[duplicates].groupby(key_cols):
                if len(group) <= 1:
                    continue
                if duplicate_handling == 'nlogp_keep_highest':
                    abs_vals = np.abs(group['nlogp'])
                    keep_idx = abs_vals.idxmax()
                    drop_idx = [idx for idx in group.index if idx != keep_idx]
                    to_drop.extend(drop_idx)
                elif duplicate_handling == 'nlogp_compare':
                    g = group.copy()
                    g['abs_nlogp'] = np.abs(g['nlogp'])
                    g_sorted = g.sort_values('abs_nlogp', ascending=False)
                    top = g_sorted.iloc[0]['abs_nlogp']
                    second = g_sorted.iloc[1]['abs_nlogp'] if len(g_sorted) > 1 else 0.0
                    rel_diff = (top - second) / top if top > 0 else 0.0
                    if rel_diff >= float(nlogp_threshold):
                        keep_idx = g_sorted.index[0]
                        drop_idx = [idx for idx in g_sorted.index[1:]]
                        to_drop.extend(drop_idx)
                    else:
                        # difference too small -> drop all in this duplicate group
                        to_drop.extend(list(g_sorted.index))
                else:
                    # Unknown strategy -> remove duplicates
                    to_drop.extend(list(group.index))

            if len(to_drop) > 0:
                before = len(go_all_merged)
                go_all_merged = go_all_merged.drop(index=list(set(to_drop)))
                after = len(go_all_merged)
                print(f"Resolved duplicates: removed {before - after} rows")
            else:
                print("No rows removed after duplicate resolution")
    else:
        print("No duplicate records found")

    # Save complete data
    go_all_merged.to_csv(output_path, index=False)

    # Save data with CC terms removed
    go_all_merged_no_cc = go_all_merged[go_all_merged['source'] != 'GO:CC'].copy()
    output_path_no_cc = output_path.replace('.csv', '_rm_CC.csv')
    go_all_merged_no_cc.to_csv(output_path_no_cc, index=False)

    print(f"Merged data saved to: {output_path}")
    print(f"Merged data (no CC) saved to: {output_path_no_cc}")
    print(f"Total records: {len(go_all_merged)}")
    print(f"Records with mapping: {go_all_merged['Group'].notna().sum()}")
    print(f"Records without mapping: {go_all_merged['Group'].isna().sum()}")

    return go_all_merged_no_cc


def predict_terms(df, classifiers, neurotransmitter_col="Neurotransmitter", go_col="GO"):
    """
    Predict GO terms for all records by subgroup (NN and N),
    while later combination still prioritizes existing mapping.

    Args:
        df: DataFrame containing GO terms
        classifiers: Dictionary of trained classifiers
        neurotransmitter_col: Neurotransmitter column name
        go_col: GO term column name

    Returns:
        DataFrame with prediction results added in column 'predicted_group'
    """
    print("=== Step 3: Predicting GO terms for all records (by Neurotransmitter) ===")

    # Copy dataframe
    df_result = df.copy()

    # Preprocess GO terms
    df_result['clean_go_term'] = df_result[go_col].apply(preprocess_go_term)

    # Initialize prediction column
    df_result['predicted_group'] = None

    if (not 'NN' in classifiers) and (not 'N' in classifiers):
        classifiers['NN'] = classifiers.get('classifier', None)
        classifiers['N'] = classifiers.get('classifier', None)
        classifiers['NN_num_to_str_mapping'] = classifiers.get('num_to_str_mapping', {})
        classifiers['N_num_to_str_mapping'] = classifiers.get('num_to_str_mapping', {})
    # Predict for NN and N groups separately on ALL rows
    if neurotransmitter_col in df_result.columns:
        # NN group prediction (all NN rows)
        nn_mask = (df_result[neurotransmitter_col] == 'NN')
        if nn_mask.any() and 'NN' in classifiers:
            print(f"Predicting {nn_mask.sum()} NN samples...")
            nn_terms = df_result.loc[nn_mask, 'clean_go_term'].values
            nn_predictions = classifiers['NN'].predict(nn_terms)

            # Convert to string labels
            if 'NN_num_to_str_mapping' in classifiers:
                nn_mapping = classifiers['NN_num_to_str_mapping']
                nn_predictions_str = [nn_mapping.get(pred, f'Unknown_{pred}') for pred in nn_predictions]
            else:
                nn_predictions_str = nn_predictions

            df_result.loc[nn_mask, 'predicted_group'] = nn_predictions_str

        # N group prediction (all N rows)
        n_mask = (df_result[neurotransmitter_col] != 'NN')
        if n_mask.any() and 'N' in classifiers:
            print(f"Predicting {n_mask.sum()} N samples...")
            n_terms = df_result.loc[n_mask, 'clean_go_term'].values
            n_predictions = classifiers['N'].predict(n_terms)

            # Convert to string labels
            if 'N_num_to_str_mapping' in classifiers:
                n_mapping = classifiers['N_num_to_str_mapping']
                n_predictions_str = [n_mapping.get(pred, f'Unknown_{pred}') for pred in n_predictions]
            else:
                n_predictions_str = n_predictions

            df_result.loc[n_mask, 'predicted_group'] = n_predictions_str

    # Handle any rows without predictions (fallback)
    no_pred_mask = df_result['predicted_group'].isna()
    if no_pred_mask.any():
        print(f"Using fallback prediction for {no_pred_mask.sum()} samples...")
        if 'NN' in classifiers:
            fallback_classifier = classifiers['NN']
            fallback_mapping = classifiers.get('NN_num_to_str_mapping', {})
        elif 'N' in classifiers:
            fallback_classifier = classifiers['N']
            fallback_mapping = classifiers.get('N_num_to_str_mapping', {})
        else:
            print("No suitable classifier found for fallback")
            fallback_classifier = None

        if fallback_classifier is not None:
            fallback_terms = df_result.loc[no_pred_mask, 'clean_go_term'].values
            fallback_predictions = fallback_classifier.predict(fallback_terms)
            fallback_predictions_str = [fallback_mapping.get(pred, f'Unknown_{pred}') for pred in fallback_predictions]
            df_result.loc[no_pred_mask, 'predicted_group'] = fallback_predictions_str

    return df_result


def evaluate_predictions(df, neurotransmitter_col="Neurotransmitter"):
    """
    Evaluate the predicted_group against the existing Group as ground truth.

    Returns a DataFrame containing overall and per-Neurotransmitter accuracy,
    as well as a per-row evaluation.
    """
    print("\n=== Evaluating predictions against mapping (ground truth) ===")

    df_eval = df.copy()

    # Only evaluate samples that have a ground-truth mapping
    mask_gt = df_eval['Group'].notna()
    num_gt = mask_gt.sum()
    if num_gt == 0:
        print("Warning: No records with ground truth found. Skipping evaluation.")
        return None, None

    df_gt = df_eval.loc[mask_gt, [
        'Region subclass' if 'Region subclass' in df_eval.columns else None,
        'Sex' if 'Sex' in df_eval.columns else None,
        'GO', 'Group', 'predicted_group', neurotransmitter_col
    ]].copy()

    # Remove placeholder None columns
    df_gt = df_gt[[c for c in df_gt.columns if c is not None]]

    # Compare predicted vs. true group
    df_gt['is_correct'] = (df_gt['predicted_group'] == df_gt['Group'])

    # Overall accuracy
    overall_acc = float(df_gt['is_correct'].mean()) if len(df_gt) > 0 else float('nan')

    # Accuracy grouped by Neurotransmitter
    if neurotransmitter_col in df_gt.columns:
        by_nt = df_gt.groupby(neurotransmitter_col)['is_correct'].agg(['mean', 'count']).reset_index()
        by_nt.rename(columns={'mean': 'accuracy', 'count': 'n_samples'}, inplace=True)
    else:
        by_nt = pd.DataFrame(columns=[neurotransmitter_col, 'accuracy', 'n_samples'])

    # Build summary DataFrame
    summary_rows = [
        {'scope': 'overall', 'subset': 'all', 'accuracy': overall_acc, 'n_samples': int(len(df_gt))}
    ]
    for _, row in by_nt.iterrows():
        summary_rows.append({
            'scope': 'by_neurotransmitter',
            'subset': row[neurotransmitter_col],
            'accuracy': float(row['accuracy']),
            'n_samples': int(row['n_samples'])
        })

    summary_df = pd.DataFrame(summary_rows)

    print(f"Ground-truth records: {len(df_gt)}")
    print("Accuracy summary:")
    print(summary_df)

    return summary_df, df_gt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GO term mapping and prediction integration script")

    # Data path parameters
    parser.add_argument("--go_dir_path", type=str,
                        default="/data2st1/junyi/output/atac1112/dar/celltype.L2/MASTNG_dar_annotated/CP_enrichment_results2/",
                        help="GO enrichment results directory path")
    parser.add_argument("--go_group_path", type=str,
                        default="/data2st1/junyi/output/atac1112/tobias/2-4 DEGl1.xlsx",
                        help="GO grouping mapping file path")
    parser.add_argument("--df_class_path", type=str,
                        default="/data2st2/junyi/code/sn/name_form_new.xlsx",
                        help="Class mapping file path")

    # Output parameters
    parser.add_argument("--output_path", type=str,
                        default="/data2st1/junyi/output/atac1112/dar/celltype.L2/MASTNG_dar_annotated/test/2_go_merged_with_predictions.csv",
                        help="Final output file path")

    # Classifier parameters
    parser.add_argument("--input_col", type=str, default='GO', help="Input feature column name (GO terms)")
    parser.add_argument("--output_col", type=str, default='Group', help="Output label column name (modules/categories)")
    parser.add_argument("--neurotransmitter_col", type=str, default="Neurotransmitter",
                        help="Neurotransmitter column name")

    # Duplicate handling parameters
    parser.add_argument("--duplicate_handling", type=str, default='filter',
                        choices=['filter', 'nlogp_compare', 'nlogp_keep_highest','Up','Down'],
                        help="How to handle duplicates in merging")
    parser.add_argument("--nlogp_threshold", type=float, default=0.1,
                        help="Relative difference threshold for abs(nlogp) when using 'nlogp_compare'")

    args = parser.parse_args()

    print("=" * 60)
    print("GO term mapping and prediction integration script")
    print("=" * 60)

    # Step 1: Merge GO data and perform mapping
    merged_df = merge_go_data(
        go_dir_path=args.go_dir_path,
        go_group_path=args.go_group_path,
        df_class_path=args.df_class_path,
        output_path=args.output_path.replace('.csv', '_merged_only.csv'),
        duplicate_handling=args.duplicate_handling,
        nlogp_threshold=args.nlogp_threshold
    )

    # Step 2: Train classifiers (using all data from go_group_path as training set)
    print("\n=== Step 2: Training classifiers ===")

    # Directly use go_group_path file to train classifiers
    classifiers = train_go_classifier(
        csv_file_path=args.go_group_path,
        input_col=args.input_col,
        output_col=args.output_col,
        neurotransmitter_col=args.neurotransmitter_col
    )

    # Step 3: Predict unmapped GO terms
    predicted_df = predict_terms(
        df=merged_df,
        classifiers=classifiers,
        neurotransmitter_col=args.neurotransmitter_col,
        go_col=args.input_col
    )

    # Step 4: Combine results
    print("\n=== Step 4: Combining mapping and prediction results ===")

    # add nlogp
    predicted_df['nlogp'] = -np.log10(predicted_df['p_value'])
    if 'Direction' in predicted_df.columns:
        predicted_df['nlogp'] = predicted_df.apply(
            lambda row: -row['nlogp'] if row['Direction'] == 'Down' else row['nlogp'],
            axis=1)

    # Create final group column: use mapping if available, otherwise use prediction
    predicted_df['final_group'] = predicted_df['Group'].fillna(predicted_df['predicted_group'])

    # Create source type column
    predicted_df['source_type'] = predicted_df['Group'].apply(lambda x: 'mapping' if pd.notna(x) else 'prediction')

    # Statistics
    mapping_count = (predicted_df['source_type'] == 'mapping').sum()
    prediction_count = (predicted_df['source_type'] == 'prediction').sum()

    print(f"Combination results:")
    print(f"- Records with mapping: {mapping_count}")
    print(f"- Records with prediction: {prediction_count}")
    print(f"- Total records: {len(predicted_df)}")

    final_df = predicted_df

    # Step 5: Save final results
    final_df.to_csv(args.output_path, index=False)
    print(f"\n=== Final results saved to: {args.output_path} ===")

    # Step 5.1: Evaluate prediction accuracy where mapping exists
    eval_summary, eval_rows = evaluate_predictions(final_df, neurotransmitter_col=args.neurotransmitter_col)
    if eval_summary is not None:
        eval_dir = os.path.dirname(args.output_path)
        eval_summary_path = args.output_path.replace('.csv', '_prediction_accuracy_summary.csv')
        eval_rows_path = args.output_path.replace('.csv', '_prediction_accuracy_rows.csv')
        eval_summary.to_csv(eval_summary_path, index=False)
        eval_rows.to_csv(eval_rows_path, index=False)
        print(f"Evaluation summary saved to: {eval_summary_path}")
        print(f"Evaluation detail rows saved to: {eval_rows_path}")

    # save to excel with two sheets: N and NN
    # Define columns to keep                Neurotransmitter
    columns_to_keep = ['Region subclass', 'Neurotransmitter', 'Class', 'Sex', 'GO', 'Group', 'predicted_group',
                       'final_group', 'source_type'
                       ]
    available_columns = [col for col in columns_to_keep if col in final_df.columns]

    df_selected = final_df[available_columns].copy()
    df_selected = df_selected.drop_duplicates()

    df_n = df_selected[df_selected['Neurotransmitter'] != 'NN'].copy()
    df_nn = df_selected[df_selected['Neurotransmitter'] == 'NN'].copy()
    output_path = args.output_path.replace('.csv', '_by_group.xlsx')
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        if len(df_n) > 0:
            df_n.to_excel(writer, sheet_name='N_group', index=False)

        if len(df_nn) > 0:
            df_nn.to_excel(writer, sheet_name='NN_group', index=False)

    # Save classifier models
    model_dir = os.path.dirname(args.output_path)
    if 'NN' in classifiers:
        nn_model_path = args.output_path.replace('.csv', '_NN_classifier.joblib')
        joblib.dump(classifiers['NN'], nn_model_path)
        print(f"NN classifier saved to: {nn_model_path}")

    if 'N' in classifiers:
        n_model_path = args.output_path.replace('.csv', '_N_classifier.joblib')
        joblib.dump(classifiers['N'], n_model_path)
        print(f"N classifier saved to: {n_model_path}")

    # Save label mappings
    mapping_data = {}
    for key in ['NN_str_to_num_mapping', 'NN_num_to_str_mapping', 'N_str_to_num_mapping', 'N_num_to_str_mapping']:
        if key in classifiers:
            mapping_data[key] = classifiers[key]

    if mapping_data:
        mapping_file = args.output_path.replace('.csv', '_label_mapping.json')
        with open(mapping_file, 'w', encoding='utf-8') as f:
            json.dump(mapping_data, f, indent=2, ensure_ascii=False)
        print(f"Label mapping saved to: {mapping_file}")

    # Print final statistics
    print(f"\n=== Final Statistics ===")
    print(f"Total records: {len(final_df)}")
    print(f"Records with original mapping: {(final_df['source_type'] == 'mapping').sum()}")
    print(f"Records with predictions: {(final_df['source_type'] == 'prediction').sum()}")
