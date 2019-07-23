#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:03:53 2019

Returns the primary anatomical structure for a gifti object.

INPUT:
G:				Nibabel gifti object

OUTPUT:
anatStruct:		AnatomicalStructurePrimary attribute from gifti object

@author: jdiedrichsen (Python conversion: switt)
"""

def getGiftiAnatomicalStruct(G):
    N = len(G._meta.data)
    anatStruct = []
    for i in range(N):
        if 'AnatomicalStructurePrimary' in G._meta.data[i].name:
            anatStruct.append(G._meta.data[i].value)
    return anatStruct
            