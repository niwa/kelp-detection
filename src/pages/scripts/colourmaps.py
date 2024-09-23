import pandas
import matplotlib.colors

def get_colourmap(dataframe: pandas.DataFrame, column_name: str, selection: str):
    """ Create a listed colourmap based on the current statuses present in the dataframe.
       Such that:
       completed = forest green
       failed = firebrick
       running = cornflowerblue
       queued / incomplete = darkturquoise
       not started = grey
    """
    colourlist = []
    if dataframe[column_name].str.contains("completed").any():
        colourlist.append("forestgreen")
    if dataframe[column_name].str.contains("failed").any():
        colourlist.append("firebrick")
    if dataframe[column_name].str.contains("incomplete").any():
        colourlist.append("magenta")
    if dataframe[column_name].str.contains("not started").any():
        colourlist.append("grey")
    if dataframe[column_name].str.contains("queued").any():
        colourlist.append("darkturquoise")
    if dataframe[column_name].str.contains("removed").any():
        colourlist.append("forestgreen")
    if dataframe[column_name].str.contains("running").any():
        colourlist.append("royalblue")
    colourlist.append("test")
    colourmap = matplotlib.colors.ListedColormap(colourlist)
    return colourmap