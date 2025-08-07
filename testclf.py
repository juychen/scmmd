# %%
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
#import gradient_boosting as gb
from sklearn.ensemble import GradientBoostingClassifier
# IMPORT support vector machines
from sklearn.svm import SVC
#from sklearn.naive_bayes import MultinomialNB
from sklearn.ensemble import RandomForestClassifier
# from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import re
import argparse
import joblib

parse = argparse.ArgumentParser(description="Train a classifier to predict module and function from GO terms")
parse.add_argument("--csv_source_path", type=str, default="/home/junyichen/code/test/GOclassifier4.csv", help="Path to the CSV file containing GO terms and their corresponding modules and functions")
parse.add_argument("--csv_test_path", type=str, default="test.csv", help="Path to the CSV file containing test GO terms")
parse.add_argument("--output_path", type=str, default=None, help="Path to save the predictions")
args = parse.parse_args()


def preprocess_go_term(go_term):
    """Clean GO term by removing prefixes and special characters"""
    # Remove GO:BP_, GO:MF_, GO:CC_ prefixes if present
    clean_term = re.sub(r'^GO:[A-Z]+_', '', go_term)
    # Remove any remaining special characters and make lowercase
    clean_term = re.sub(r'[^\w\s]', '', clean_term).lower()
    return clean_term

def train_go_classifier(csv_file_path):
    """
    Train a classifier to predict module and function from GO terms
    
    Args:
        csv_file_path: Path to CSV file with columns: go_term, module, function
        
    Returns:
        A trained pipeline that can predict module and function for new GO terms
    """
    # Load and preprocess data
    df = pd.read_csv(csv_file_path)
    df['clean_go_term'] = df.iloc[:, 0].apply(preprocess_go_term)
    # Find nan from df
    df = df.dropna()
    # downsample  class 10 in df to
    # Split into features (X) and targets (y)
    X = df['clean_go_term']
    y_module = df.iloc[:, 1]  # Module number
    y_function = df.iloc[:, 2]  # Function
    
    # Split data into training and test sets
    X_train, X_test, y_module_train, y_module_test, y_function_train, y_function_test = train_test_split(
        X, y_module, y_function, test_size=0.2, random_state=42
    )
    
    # Create pipeline for module classification
    # module_pipeline = Pipeline([
    #     ('tfidf', TfidfVectorizer(ngram_range=(1, 2), max_features=5000)),
    #     ('clf', LogisticRegression(multi_class='multinomial', solver='lbfgs', max_iter=1000))
    # ])
    
    # Create pipeline for function classification
    module_pipeline = Pipeline([
        ('tfidf', TfidfVectorizer(ngram_range=(1, 2), max_features=5000)),
        #('clf', GradientBoostingClassifier(n_estimators=100, random_state=42))  # Using Gradient Boosting for module classification
        #('clf', RandomForestClassifier(n_estimators=100, random_state=42))  # Using Random Forest for module classification
        ('clf', SVC(kernel='linear'))  # Using SVM for function classification
    ])
    
    # Train both classifiers
    module_pipeline.fit(X_train, y_module_train)
    # module_svc.fit(X_train, y_module_train)
    
    # Evaluate performance
    module_pred = module_pipeline.predict(X_test)
    # svc_pred = module_svc.predict(X_test)
    #function_pred = function_pipeline.predict(X_test)
    
    print("Module Classification Report:")
    print(classification_report(y_module_test, module_pred))
    
    # print("\nFunction Classification Report:")
    # print(classification_report(y_module_test, svc_pred))
    
    # Return trained models
    return {
        'module_classifier': module_pipeline
        # 'function_classifier': module_svc
    }

def predict_go_term(classifiers, go_term):
    """
    Predict module and function for a new GO term
    
    Args:
        classifiers: Dictionary containing 'module_classifier' and 'function_classifier'
        go_term: New GO term to classify
        
    Returns:
        Dictionary with predicted module and function
    """
    # Preprocess the input GO term
    clean_term = preprocess_go_term(go_term)
    
    # Make predictions
    module = classifiers['module_classifier'].predict([clean_term])[0]
    # function = classifiers['function_classifier'].predict([clean_term])[0]
    
    return {
        'go_term': go_term,
        'predicted_module': module
        # 'predicted_function': function
    }


# Replace with your actual CSV file path
csv_path = args.csv_source_path
csv_test = args.csv_test_path
# Train the classifiers
classifiers = train_go_classifier(csv_path)

# # Save the trained models for later use
joblib.dump(classifiers['module_classifier'], 'go_module_classifier.joblib')
# joblib.dump(classifiers['function_classifier'], 'go_function_classifier.joblib')

# # Example prediction
# new_go_term = "GO:0005740	GO:CC	mitochondrial envelope"
# prediction = predict_go_term(classifiers, new_go_term)
# print("\nPrediction for new GO term:")
# print(prediction)


# Load the trained classifiers
predictors = classifiers['module_classifier']
# joblib.dump(predictors, 'go_module_classifier.joblib')
# read the test data
df_test = pd.read_csv(csv_test)
df_test['clean_go_term'] = df_test.iloc[:, 0].apply(preprocess_go_term)
# Predict on the test data
predictions = df_test['clean_go_term'].apply(lambda x: predict_go_term({'module_classifier': predictors}, x))   
# Convert predictions to DataFrame
predictions_df = pd.DataFrame(predictions.tolist())
# Save predictions to CSV
dfresult = pd.DataFrame(
    {'go_term': df_test.iloc[:, 0], 'predicted_module': predictions_df['predicted_module']}
)
df_test['predicted_module'] = predictions_df['predicted_module']
if args.output_path is None:
    df_test.to_csv(args.csv_test_path+".predicted.csv", index=False)
else:
    df_test.to_csv(args.output_path, index=False)

