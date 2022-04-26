# All models from here:
#http://prd.stsci.edu/prd/files/UVMelem.cgi?ELEM=saa-models%2Fmodel31.svdf&DB=ops&REV=1.4
# Starting with 1.1, through 1.4
# Tuples are latitude, longitude

#**** SAA MODEL 31 COS FUV (XDL)                                        GC
#     INITIAL ENTRY: CLONE OF MODEL 25.  PR 45488,       03/22/02       GC
model_0 = [ 
(-28.3,  20.0),
(-27.5,  21.0),
(-26.1,  19.0),
(-19.8,   7.5),
( -9.6, 347.0),
( -7.6, 336.4),
( -6.0, 324.8),
( -7.9, 303.2),
(-12.0, 292.1),
(-17.1, 285.9),
(-20.3, 283.5),
(-23.5, 282.5),
(-26.0, 282.4),
(-28.6, 282.7),
(-28.3,  20.0)]

#**** SAA MODEL 31 COS FUV (XDL)                                        GC
#     INITIAL ENTRY: CLONE OF MODEL 25.  PR 45488,       03/22/02       GC
#     UPDATE:  SHIFTED WEST 6DEG.  PR 65147, 5/13/10                    GC
model_1 = [ 
(-28.3,  14.0),
(-27.5,  15.0),
(-26.1,  13.0),
(-19.8,   1.5),
( -9.6, 341.0),
( -7.6, 330.4),
( -6.0, 318.8),
( -7.9, 297.2),
(-12.0, 286.1),
(-17.1, 279.9),
(-20.3, 277.5),
(-23.5, 276.5),
(-26.0, 276.4),
(-28.6, 276.7),
(-28.3,  14.0)]

#**** SAA MODEL 31 COS FUV (XDL)                                        GC
#     INITIAL ENTRY: CLONE OF MODEL 25.  PR 45488,       03/22/02       GC
#     UPDATE:  SHIFTED WEST 6DEG.  PR 65147, 5/13/10                    GC
#     UPDATE: SHIFT NW and N.  PR 69747, 11/8/11                        GC
model_2 = [ 
(-28.3,  14.0),
(-27.5,  15.0),
(-26.1,  13.0),
(-19.8,   1.5),
( -9.6, 341.0),
( -5.5, 330.4),
( -3.5, 318.8),
( -3.8, 308.0),
( -6.2, 295.5),
(-10.5, 285.6),
(-16.0, 279.5),
(-20.3, 277.3),
(-23.5, 276.5),
(-26.0, 276.4),
(-28.6, 276.7),
(-28.3,  14.0)]

#**** SAA MODEL 31 COS FUV (XDL)                                        GC
#     INITIAL ENTRY: CLONE OF MODEL 25.  PR 45488,       03/22/02       GC
#     UPDATE:  SHIFTED WEST 6DEG.  PR 65147, 5/13/10                    GC
#     UPDATE: SHIFT NW and N.  PR 69747, 11/8/11                        GC
#     Update: NW expansion.  PDB-91385 9/6/19                           GC
model_3 = [ 
(-29.0,  14.0), 
(-27.5,  15.0),
(-26.1,  13.0),
(-19.8,   1.5),
( -9.6, 341.0),
( -5.5, 330.4),
( -3.4, 318.8),
( -3.6, 308.0),
( -5.3, 295.4),
( -9.6, 284.8),
(-15.0, 278.4),
(-20.3, 276.4),
(-23.5, 276.1),
(-26.0, 276.0),
(-29.0, 276.4),
(-29.0,  14.0)] 

saa_models = {(47892, 55329): model_0,
              (55329, 55873): model_1,
              (55873, 58732): model_2,
              (58732, 60676): model_3}

def get_saa_poly(mjd):
    for drange in saa_models:
        if mjd > drange[0] and mjd <= drange[1]:
            model = saa_models[drange]
    saa_lat = [x[0] for x in model]
    saa_lon = [x[1] for x in model]
    return saa_lat, saa_lon

