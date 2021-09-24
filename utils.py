import pandas as pd

def sql_to_df(sql_results, returncols):
    # Convert query results (list of attributes) to pandas dataframe
#    d = {col: [getattr(x, col) for x in results] for col in returncols}
    d = {}
    for col in returncols:
        try:
            d[col] = [getattr(x, col) for x in sql_results]
        except:
            pass
    df = pd.DataFrame(d)

    return df

def timefunc(func):  # IN USE
    """
    Decorator to wrap and time functions.

    Parameters
    ----------
    func: func
        Function to be timed

    Returns
    -------
    wrapper: func
        Timed version of the function
    """

    def wrapper(*args, **kw):
        t1 = datetime.datetime.now()
        result = func(*args, **kw)
        t2 = datetime.datetime.now()

        print(f"{func.__name__} executed in {t2 - t1}")

        return result

    return wrapper


