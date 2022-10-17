import os
import sqlite3

dbt = {int:"INTEGER", float:"REAL", str:"TEXT", type(None):"NULL", bool:"INTEGER"}

def add_to_thermo_table(name, thermo_data, replace=True):
    direc = os.path.dirname(os.path.abspath(__file__))
    database_file = os.path.join(direc, "models", "testing.db")
    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(""" SELECT count(name) FROM sqlite_master WHERE type='table' AND
                name='THERMO_DATA' """)
    # Creating table if it does not exist
    if cursor.fetchone()[0] != 1:
        table = """CREATE TABLE THERMO_DATA (id INTEGER PRIMARY KEY, Name TEXT);"""
        cursor.execute(table)
    # delete old entry
    if replace:
        cursor.execute(f"DELETE FROM THERMO_DATA WHERE Name = '{name}'")
    cursor.execute(f""" SELECT count(name) FROM THERMO_DATA WHERE Name='{name}' """)
    # entry already exists
    if cursor.fetchone()[0] != 1:
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

def append_runtime_table(name, runtime):
    direc = os.path.dirname(os.path.abspath(__file__))
    database_file = os.path.join(direc, "models", "testing.db")
    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(""" SELECT count(name) FROM sqlite_master WHERE type='table' AND
                name='PERFORMANCE' """)
    # Creating table if it does not exist
    if cursor.fetchone()[0] != 1:
        table = """CREATE TABLE PERFORMANCE (id TEXT PRIMARY KEY, time REAL, nruns INTEGER);"""
        cursor.execute(table)
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

def get_database_connection(name):
    direc = os.path.dirname(os.path.abspath(__file__))
    database_file = os.path.join(direc, "models", "testing.db")
    connection = sqlite3.connect(database_file)
    return connection

def create_simulation_table(name):
    direc = os.path.dirname(os.path.abspath(__file__))
    database_file = os.path.join(direc, "models", "testing.db")
    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(f"""DROP TABLE IF EXISTS {name}""")
    cursor.execute(f"""CREATE TABLE {name} (id INTEGER PRIMARY KEY);""")
    connection.commit()

def first_time_step(connection, name, stats):
    cursor = connection.cursor()
    cursor.execute(f"""INSERT INTO {name} ( id ) VALUES ( 1 )""")
    for k, v in stats.items():
        try:
            cursor.execute(f"ALTER TABLE {name} ADD COLUMN {k} {dbt[type(v)]}")
        except Exception as e:
            pass
    for k, v in stats.items():
            value = f"\'{v}\'" if isinstance(v, str) else v
            cursor.execute(f"""UPDATE {name} SET {k} = {value} WHERE id = 1 """)
    connection.commit()
    return 1

def add_time_step(connection, id, name, stats):
    cursor = connection.cursor()
    cursor.execute(f"""INSERT INTO {name} ( id ) VALUES ( {id} )""")
    for k, v in stats.items():
            value = f"\'{v}\'" if isinstance(v, str) else v
            cursor.execute(f"""UPDATE {name} SET {k} = {value} WHERE id = {id} """)
    connection.commit()
