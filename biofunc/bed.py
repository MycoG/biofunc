import pandas as pd

def load_bed(path:str, header=True) -> pd.DataFrame:
    """
    load BED file as pandas DataFrame

    :param path: path to BED file
    :param header: look for .header file in same directory
    :type header: bool
    """
    #TODO: allow for headers to be added, like pd.read_csv(names=)

    if header:
        with open(path+".header", 'r') as f:
            cols = f.read().strip().split("\n")
        df = pd.read_csv(path, names=cols, sep='\t', index_col=False)
        return df
    
    else:
        df = pd.read_csv(path, sep='\t', index_col=False, header=None)
        return df

def save_bed(df:pd.DataFrame, out:str) -> None:
    """
    save pandas dataframe as a BED file and header file

    :param df: dataframe
    :type df: pd.DataFrame
    :param out: output filename
    :type out: str
    """
    df.to_csv(out, sep="\t", index=False, header=False)
    with open(out+".header", 'w') as f:
        f.write("\n".join(list(df.columns)) + "\n")
    return