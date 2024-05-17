import pandas as pd


if __name__ == "__main__":
    predictions = pd.read_csv("/home/malte/projects/refactorring/majority_test/main.tsv",
                              sep="\\t", header=0, names=['mirna', 'target', 'start', 'end', 'tool' ])

    start = False
    end = True
    complete = False  # TODO: Ask if this means an exact match or an "either ..., or ..."
    majority = 2

    if start:  # group by start indices
        predictions = predictions.groupby(['mirna', 'target', 'start'])['tool'].apply(set).reset_index()
    elif end:  # group by end indices
        predictions = predictions.groupby(['mirna', 'target', 'end'])['tool'].apply(set).reset_index()
    elif complete:  # group by both indices
        predictions = predictions.groupby(['mirna', 'target', 'start', 'end'])['tool'].apply(set).reset_index()

    # performing majority vote keeping only mirna binding sites that meet the required number of votes
    post_vote_predictions = predictions[predictions['tool'].apply(len) >= majority].copy()
    post_vote_predictions = post_vote_predictions.drop('tool', axis=1)


    predictions_out.to_csv('filtered_bindingsites.tsv', sep='\\t', index=False)
