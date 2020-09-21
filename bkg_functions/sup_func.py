# -*- coding: utf-8 -*-
"""
Created on Sunday September 20 12:23:08 2020
@author: pcaldas
"""
import pandas as pd
import numpy as np
from xml.etree import cElementTree as ET

def read_trackmate_xml_tracks(xml_file):
    """Reads tracks from trackmate xml track file 
    and converts into a user-friendly DataFrame """
	
    tracks = ET.parse(xml_file)

    frame_interval = float(tracks.getroot().attrib["frameInterval"])

    n_tracks = float(tracks.getroot().attrib["nTracks"])

    attributes = []

    for ti, track in enumerate(tracks.iterfind('particle')):
        for spots in track.iterfind('detection'):
            attributes.append([ti, int(spots.attrib.get('t')),
                           float(spots.attrib.get('x')),
                           float(spots.attrib.get('y'))])

    track_table = pd.DataFrame(attributes, columns=['TRACK_ID','FRAME','POSITION_X','POSITION_Y'])

    return track_table, frame_interval, n_tracks