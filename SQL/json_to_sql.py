import pandas as pd
import mysql.connector
from sqlalchemy import create_engine
from sqlalchemy import text
import json
import argparse

# mydb = mysql.connector.connect(
#     host="localhost",
#     port=3307,
#     user="root",
#     password="pixhawk2",
#     database="cas_delay_test"
#     )

engine = create_engine("mysql+mysqlconnector://root:pixhawk2@localhost:3307/cas_delay_test")


def to_sql(spacing, exptype, scen):

    if 'Bez' in exptype or 'Bez' in scen:
        '''
        Bez Data
        '''
        file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Bez_{spacing}_Apart.json"
        with open(file_path, 'r') as file:
            BezierData = json.load(file)
        records = []
        for item in BezierData:
            for key, value in item.items():
                records.append({
                    # 'ID': key,
                    'p0x': value['p0x'],
                    'p0y': value['p0y'],
                    'p1x': value['p1x'],
                    'p1y': value['p1y'],
                    'p2x': value['p2x'],
                    'p2y': value['p2y'],
                    'length': value['length'],
                    'Bez_ID': value['Bez_ID'],
                    'TOA': value['TOA'],
                    'orig_TOA': value['orig_TOA'],
                    'travel_time': value['travel_time'],
                    'timeStamp': value['timeStamp'],
                    'velocity': value['velocity'],
                    'ACID': value['ACID'],
                    'ExpNum': value['ExpNum'],
                    'Category': value['Category'],
                    'ExpType': value['ExpType'],
                    'ax_spacing': value['ax_spacing']
                })

        bez = pd.DataFrame(records)
        bez.to_sql(name = 'bez', con=engine, if_exists='append', index = False)
        print('Bez Data Uploaded to SQL!')


        '''
        Dubins Data

        '''
        file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Dubins_{spacing}_Apart.json"

        with open(file_path, 'r') as file:
            DubinsData = json.load(file)

        # Process JSON into records
        records = []
        for item in DubinsData:
            for path_type, value in item.items():
                records.append({
                    'path_type': value.get('path_type'),
                    'bez_intx': value.get('bez_intx'),
                    'bez_inty': value.get('bez_inty'),
                    'nom_intx': value.get('nom_intx'),
                    'nom_inty': value.get('nom_inty'),
                    'bez_t': value.get('bez_t'),
                    'intersect_heading': value.get('intersect_heading'),
                    'h': value.get('h'),
                    'k': value.get('k'),
                    'bank_angle': value.get('bank_angle'),
                    'tr': value.get('tr'),
                    'timeStamp': value.get('timeStamp'),
                    'ACID': value.get('ACID'),
                    'ExpNum': value.get('ExpNum'),
                    'Category': value.get('Category'),
                    'ExpType': value.get('ExpType'),
                    'ax_spacing': value.get('ax_spacing'),
                })


        dubinsPath = pd.DataFrame(records)
        dubinsPath.to_sql(name='dubins', con=engine, if_exists='append', index=False)
        print('Dubins Data Uploaded to SQL!')

    '''
    State Data
    '''  
    file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\State_{spacing}_Apart.json"

    with open(file_path) as file:
        state_data = json.load(file)

    records = []
    for key, value in state_data.items():
        for i in range(len(value['timeStamp'])):
            records.append({
                # 'ID': key,
                'ACID': value['ACID'],
                'ExpNum': value['ExpNum'],
                'Category': value['Category'],
                'ExpType': value['ExpType'],
                'timeStamp': value['timeStamp'][i],
                'lat': value['lat'][0][i],
                'lon': value['lon'][0][i],
                'Heading': value['Heading'][0][i],
                'V_meters': value['V_meters'][0][i],
                'V_knots': value['V_knots'][0][i],
                'ax_spacing': value['ax_spacing']
            })

    state = pd.DataFrame(records)
    state.to_sql(name = 'state', con=engine, if_exists='append', index = False)
    print('State Data Uploaded to SQL!')

    '''
    EV Specific
    '''
    file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\EVSpecific_{spacing}_Apart.json"
    with open(file_path, 'r') as file:
        ev_data = json.load(file)

    # Convert the dictionary to a pandas DataFrame
    records = []
    for value in ev_data:
        records.append({
            'TOI': value['TOI'],  
            'TOI_Dist': value['TOI_Dist'],  
            'timeStamp': value['timeStamp'],  
            'ax_spacing': value['ax_spacing'],  
            'ACID': value['ACID'],  
            'ExpNum': value['ExpNum'],  
            'Category': value['Category'],  
            'ExpType': value['ExpType']
        })

    # Create DataFrame
    ev_specific = pd.DataFrame(records)
    ev_specific.to_sql(name = 'ev_specific', con=engine, if_exists='append', index = False)
    print('EV Data Uploaded to SQL!')

    '''
    Notable Events
    '''
    file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\NotableEvents_{spacing}_Apart.json"
    with open(file_path, 'r') as file:
        note_events = json.load(file)

    # Convert the dictionary to a pandas DataFrame
    record2 = []
    for value in note_events:
        record2.append({
            'event': value['event'],  
            'timeStamp': value['timeStamp'],  
            'ACID': value['ACID'],  
            'ExpNum': value['ExpNum'],  
            'Category': value['Category'],  
            'ExpType': value['ExpType'],
            'ax_spacing': value['ax_spacing']
        })

    # Create DataFrame
    notable_events = pd.DataFrame(record2)
    notable_events.to_sql(name = 'note_events', con=engine, if_exists='append', index = False)
    print('Event Data Uploaded to SQL!')
    # print(notable_events)

    '''
    Aircraft Data
    '''
    file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Aircraft_{spacing}_Apart.json"
    with open(file_path, 'r') as file:
        ac_dat = json.load(file)

    # Convert the dictionary to a pandas DataFrame
    records = []
    for value in ac_dat:
        records.append({
            'ACID': value['ACID'],
            'ExpNum': value['ExpNum'],
            'Category': value['Category'],
            'dt': value['dt'],
            'ExpType': value['ExpType'],
            'ax_spacing': value['ax_spacing']
        })

    # Create DataFrame
    aircraft_data = pd.DataFrame(records)
    aircraft_data.to_sql(name = 'aircraft', con=engine, if_exists='append', index = False)
    print('Aircraft Data Uploaded to SQL!')

    '''
    Delay Data
    '''
    file_path = f"C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Delay_{spacing}_Apart.json"
    with open(file_path, 'r') as file:
        DelayData = json.load(file)
    records = []

    for key, value in DelayData.items():
        records.append({
            # 'ID': key,
            'delay': value['delay'],
            'nom_time': value['nom_time'],
            'am_time': value['am_time'],
            'start_lat': value['start_lat'],
            'start_lon': value['start_lon'],
            'end_lat': value['end_lat'],
            'end_lon': value['end_lon'],
            'point_dist': value['point_dist'],
            'ACID': value['ACID'],
            'ExpNum': value['ExpNum'],
            'Category': value['Category'],
            'ExpType': value['ExpType'],
            'ax_spacing': value['ax_spacing']
        })

    delay = pd.DataFrame(records)
    delay.to_sql(name = 'delay', con=engine, if_exists='append', index = False)
    print('Delay Data Uploaded to SQL!')

if __name__ == '__main__':
    windows = True
    if not windows:
        parser = argparse.ArgumentParser(description='Apply different spacing between fleet aircraft')
        parser.add_argument('-s1', '--scenario')
        # parser.add_argument('-s2', '--subscenario')
        parser.add_argument('-d', '--spacing')
        # parser.add_argument('-t', '--time')
        # parser.add_argument('-en', '-expnum')
        parser.add_argument('-et', '--exptype')
        args = parser.parse_args()
        to_sql(args.spacing, args.exptype, args.scenario)
    else:
        to_sql( 311, "SuperTest", 'BezAM')