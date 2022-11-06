import sqlite3
import matplotlib.pyplot as plt
from cantera_adaptive_testing.database_utils import *

# use appropriate backend
# change font
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'serif'
# getting yaml for use in functions

def plot_contour_threshold_time(database, model_name, problem, prop, **kwargs):
    # create record names
    names = []
    classifiers = []
    if kwargs.pop("remove_thirdbody", False):
        classifiers.append("ntb")
    if kwargs.pop("remove_falloff", False):
        classifiers.append("nfo")
    classifiers.append(problem)
    for i in range(0, 19):
        names.append("_".join([model_name, str(i)]+classifiers))
    # get database of interest
    connect = get_database_connection(database)
    cursor = connect.cursor()
    cursor.execute("SELECT * FROM dbname.sqlite_master WHERE type='table'")
    names = cursor.fetchall()
    print(names)
    # query each table for the condition of interest
    # for n in names:
    #     query = f"SELECT {prop} FROM "+"{:s}"
    #     cursor.execute(query.format(n))
    #     res = cursor.fetchall()
    #     print(res)
    #     input()


if __name__ == "__main__":
    plot_contour_threshold_time("surf-perf-temp.db", "PlatinumSmallHydrogen", "plug_flow_reactor", "condition")
