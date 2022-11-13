import os
import re
import sqlite3
import numpy as np


dbt = {int:"INTEGER", float:"REAL", np.float64:"REAL", str:"TEXT", type(None):"NULL", bool:"INTEGER"}


def add_to_thermo_table(name, thermo_data, database):
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    # delete old entry
    cursor.execute(f""" SELECT count(name) FROM THERMO_DATA WHERE Name='{name}' """)
    if cursor.fetchone()[0] != 0:
        cursor.execute(f"DELETE FROM THERMO_DATA WHERE Name = '{name}'")
    # write new entry
    cursor.execute(f"""INSERT INTO THERMO_DATA ( Name ) VALUES ( '{name}' )""")
    for k, v in thermo_data.items():
        try:
            cursor.execute(f"ALTER TABLE THERMO_DATA ADD COLUMN {k} {dbt[type(v)]}")
        except Exception as e:
            pass
    for k, v in thermo_data.items():
        value = f"\'{v}\'" if isinstance(v, str) else v
        cursor.execute(f"""UPDATE THERMO_DATA SET {k} = {value} WHERE Name = '{name}' """)
    connection.commit()
    cursor.execute(f"SELECT id FROM THERMO_DATA WHERE Name = '{name}'")
    return cursor.fetchone()[0]


def create_all_tables(database):
    # create other tables
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    # create performance table if it doesn't exist
    table = """CREATE TABLE IF NOT EXISTS PERFORMANCE (id TEXT PRIMARY KEY, time REAL, nruns INTEGER);"""
    cursor.execute(table)
    # create exceptions table if it doesn't exist
    table = """CREATE TABLE IF NOT EXISTS EXCEPTIONS (id TEXT PRIMARY KEY, exception TEXT);"""
    cursor.execute(table)
    # create thermo_data table if it does not exist
    table = """CREATE TABLE IF NOT EXISTS THERMO_DATA (id INTEGER PRIMARY KEY, Name TEXT);"""
    cursor.execute(table)
    connection.commit()
    connection.close()


def get_steadystate_time(name):
    # see if performance data already exists
    # create other tables
    direc = os.path.dirname(os.path.abspath(__file__))
    database = os.path.join(direc, "models", "steady_state_times.db")
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    cursor.execute(f""" SELECT steadytime FROM STEADY_STATE_TIME WHERE id='{name}' """)
    return cursor.fetchone()


def append_steadystate_time_table(name, sstime):
    direc = os.path.dirname(os.path.abspath(__file__))
    database = os.path.join(direc, "models", "steady_state_times.db")
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    # create table if it doesn't exist
    table = """CREATE TABLE IF NOT EXISTS STEADY_STATE_TIME (id TEXT PRIMARY KEY, steadytime REAL);"""
    cursor.execute(table)
    # see if performance data already exists
    cursor.execute(f""" SELECT count(id) FROM STEADY_STATE_TIME WHERE id='{name}' """)
    if cursor.fetchone()[0] != 0:
        cursor.execute(f"DELETE FROM STEADY_STATE_TIME WHERE id='{name}'")
    cursor.execute(f"""INSERT INTO STEADY_STATE_TIME ( id, steadytime ) VALUES ( '{name}', {sstime} )""")
    connection.commit()


def append_runtime_table(name, runtime, database):
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    # see if performance data already exists
    cursor.execute(f""" SELECT count(id) FROM PERFORMANCE WHERE id='{name}' """)
    if cursor.fetchone()[0] == 0:
        cursor.execute(f"""INSERT INTO PERFORMANCE ( id, time, nruns) VALUES ( '{name}', {runtime}, 1 )""")
    else:
        cursor.execute(f""" SELECT * FROM PERFORMANCE WHERE id='{name}' """)
        entry = cursor.fetchone()
        new_rt = (runtime + entry[1])/2
        new_runs = entry[2] + 1
        cursor.execute(f"""UPDATE PERFORMANCE SET time = {new_rt}, nruns = {new_runs} WHERE id = '{name}'""")
    connection.commit()


def append_exception_table(name, exception, database):
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    # see if performance data already exists
    cursor.execute(f""" SELECT count(id) FROM EXCEPTIONS WHERE id='{name}' """)
    if cursor.fetchone()[0] == 0:
        cursor.execute(f"""INSERT INTO EXCEPTIONS (id, exception) VALUES ( '{name}', '{exception}' )""")
    else:
        cursor.execute(f"""UPDATE EXCEPTIONS SET exception = '{exception}' WHERE id = '{name}'""")
    connection.commit()


def create_simulation_table(name, database):
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    cursor.execute(f"""DROP TABLE IF EXISTS {name}""")
    cursor.execute(f"""CREATE TABLE {name} (id INTEGER PRIMARY KEY);""")
    connection.commit()


def add_analysis_stats(name, arr, keys, database):
    connection = sqlite3.connect(database, timeout=100000)
    cursor = connection.cursor()
    # create table first
    cursor.execute(f"""DROP TABLE IF EXISTS {name}""")
    ct_keys = ", ".join([ f"{k} {dbt[type(v)]}" for k, v in zip(keys[1:], arr[0, 1:])])
    cursor.execute(f"""CREATE TABLE {name} (id REAL PRIMARY KEY, {ct_keys});""")
    # insert rows into table
    key_str = ", ".join(keys)
    placeholders = ', '.join('?' * len(keys))
    sql = f"INSERT INTO {name} ({key_str}) VALUES ({placeholders})"
    for row in arr:
        cursor.execute(sql, row)
    connection.commit()


def get_performance_data(connection, thread):
    pass
