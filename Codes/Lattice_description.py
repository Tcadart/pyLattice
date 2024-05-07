#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Lattice geometry
#Sequence of tuples defining each member in the unit cell, with each tuple formatted (startX, startY, startZ, endX, endY, endZ, radius)
#Warning : Beam need to be define without center node
#*******************************************************************************************************************
#*******************************************************************************************************************
# 0 => BCC
# 1 => Octet
# 2 => OctetExt
# 3 => OctetInt
# 4 => BCCZ
# 5 => Cubic
# 6 => OctahedronZ
# 7 => OctahedronZcross
# 8 => Kelvin
# 9 => Cubic formulation 2 (centered)
# 10 => Cubic V3
# 11 => Cubic V4
# 12 => New lattice (non connu) GPT generated
# 13 => Diamond
def Lattice_geometry(Lattice,Radius_geom):
    BCC = [(0.0, 0.0, 0.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
        (0.5, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
        (0.5, 0.5, 0.5, 0.0, 0.0, 1.0, Radius_geom),
        (0.5, 0.5, 0.5, 0.0, 1.0, 0.0, Radius_geom),
        (0.5, 0.5, 0.5, 0.0, 1.0, 1.0, Radius_geom),
        (1.0, 0.0, 1.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 1.0, 0.0, 0.0, Radius_geom)]
    Octet = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
        (0.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
        (0.0, 0.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 1.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 0.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 1.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (1.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (1.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (0.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (0.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (0.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
        (1.0, 0.0, 0.0, 1.0, 0.5, 0.5, Radius_geom),
        (1.0, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
        (1.0, 0.0, 1.0, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 1.0, 1.0, Radius_geom),
        (0.0, 1.0, 0.0, 0.5, 1.0, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 1.0, 0.0, Radius_geom),
        (0.0, 1.0, 1.0, 0.5, 1.0, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom)]
    OctetExt = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
        (0.0, 0.0, 1.0, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 0.0, 0.5, 0.0, 0.5, Radius_geom),
        (0.0, 0.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 1.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 0.0, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 1.0, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (1.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (1.0, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (0.0, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
        (0.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (0.0, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
        (1.0, 0.0, 0.0, 1.0, 0.5, 0.5, Radius_geom),
        (1.0, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
        (1.0, 0.0, 1.0, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 1.0, 1.0, Radius_geom),
        (0.0, 1.0, 0.0, 0.5, 1.0, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 1.0, 0.0, Radius_geom),
        (0.0, 1.0, 1.0, 0.5, 1.0, 0.5, Radius_geom)]
    OctetInt = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom)]
    BCCZ = [(0.5, 0.5, 0.5, 1.0, 1.0, 1.0, Radius_geom),
        (0.0, 0.0, 0.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 1.0, 1.0, 0.0, Radius_geom),
        (0.0, 0.0, 1.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 0.0, 1.0, 0.0, Radius_geom),
        (1.0, 0.0, 1.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 0.0, 1.0, 1.0, Radius_geom),
        (1.0, 0.0, 0.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 0.5, 0.5, 1.0, Radius_geom)]
    Cubic = [(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, Radius_geom),
          (1.0, 0.0, 0.0, 1.0, 0.0, 1.0, Radius_geom),
          (0.0, 1.0, 0.0, 0.0, 1.0, 1.0, Radius_geom),
          (1.0, 1.0, 0.0, 1.0, 1.0, 1.0, Radius_geom),
          (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, Radius_geom),
          (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, Radius_geom),
          (1.0, 1.0, 0.0, 0.0, 1.0, 0.0, Radius_geom),
          (1.0, 1.0, 0.0, 1.0, 0.0, 0.0, Radius_geom),
          (0.0, 0.0, 1.0, 1.0, 0.0, 1.0, Radius_geom),
          (0.0, 0.0, 1.0, 0.0, 1.0, 1.0, Radius_geom),
          (1.0, 1.0, 1.0, 0.0, 1.0, 1.0, Radius_geom),
          (1.0, 1.0, 1.0, 1.0, 0.0, 1.0, Radius_geom)]
    OctahedronZ = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.0, 0.5, 0.5, 1.0, Radius_geom)]
    OctahedronZcross = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 1.0, 0.5, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 1.0, Radius_geom),
        (1.0, 0.5, 0.5, 0.5, 0.5, 0.0, Radius_geom),
        (0.5, 0.5, 0.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 0.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 1.0, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.0, 0.5, 0.5, 0.5, 0.5, Radius_geom),
        (0.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
        (1.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 1.0, 0.5, 0.5, 0.5, 0.5, Radius_geom)]
    Kelvin = [(0.5, 0.25, 0, 0.25, 0.5, 0, Radius_geom),
              (0.5, 0.25, 0, 0.75, 0.5, 0, Radius_geom),
              (0.5, 0.75, 0, 0.25, 0.5, 0, Radius_geom),
              (0.5, 0.75, 0, 0.75, 0.5, 0, Radius_geom),
              (0.5, 0.25, 1, 0.25, 0.5, 1, Radius_geom),
              (0.5, 0.25, 1, 0.75, 0.5, 1, Radius_geom),
              (0.5, 0.75, 1, 0.25, 0.5, 1, Radius_geom),
              (0.5, 0.75, 1, 0.75, 0.5, 1, Radius_geom),
              (0.5, 0, 0.25, 0.25, 0, 0.5, Radius_geom),
              (0.5, 0, 0.25, 0.75, 0, 0.5, Radius_geom),
              (0.5, 0, 0.75, 0.25, 0, 0.5, Radius_geom),
              (0.5, 0, 0.75, 0.75, 0, 0.5, Radius_geom),
              (0.5, 1, 0.25, 0.25, 1, 0.5, Radius_geom),
              (0.5, 1, 0.25, 0.75, 1, 0.5, Radius_geom),
              (0.5, 1, 0.75, 0.25, 1, 0.5, Radius_geom),
              (0.5, 1, 0.75, 0.75, 1, 0.5, Radius_geom),
              (0, 0.5, 0.25, 0, 0.25, 0.5, Radius_geom),
              (0, 0.5, 0.25, 0, 0.75, 0.5, Radius_geom),
              (0, 0.5, 0.75, 0, 0.25, 0.5, Radius_geom),
              (0, 0.5, 0.75, 0, 0.75, 0.5, Radius_geom),
              (1, 0.5, 0.25, 1, 0.25, 0.5, Radius_geom),
              (1, 0.5, 0.25, 1, 0.75, 0.5, Radius_geom),
              (1, 0.5, 0.75, 1, 0.25, 0.5, Radius_geom),
              (1, 0.5, 0.75, 1, 0.75, 0.5, Radius_geom),
              (0.5, 0.25, 0, 0.5, 0, 0.25, Radius_geom),
              (0.25, 0.5, 0, 0, 0.5, 0.25, Radius_geom),
              (0.75, 0.5, 0, 1, 0.5, 0.25, Radius_geom),
              (0.5, 0.75, 0, 0.5, 1, 0.25, Radius_geom),
              (0.25, 0, 0.5, 0, 0.25, 0.5, Radius_geom),
              (0.75, 0, 0.5, 1, 0.25, 0.5, Radius_geom),
              (0.75, 1, 0.5, 1, 0.75, 0.5, Radius_geom),
              (0.25, 1, 0.5, 0, 0.75, 0.5, Radius_geom),
              (0.5, 0, 0.75, 0.5, 0.25, 1, Radius_geom),
              (0, 0.5, 0.75, 0.25, 0.5, 1, Radius_geom),
              (0.5, 1, 0.75, 0.5, 0.75, 1, Radius_geom),
              (1, 0.5, 0.75, 0.75, 0.5, 1, Radius_geom)]
    CubicV2 = [(0.5, 0.0, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (0.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 1.0, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (1.0, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 0.5, 0.0, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 0.5, 1.0, 0.5, 0.5, 0.5, Radius_geom)]
    CubicV3 = [(0.0, 0.0, 0.0, 0.5, 0.0, 0.0, Radius_geom),
                (0.5, 0.0, 0.0, 1.0, 0.0, 0.0, Radius_geom),
                (0.0, 1.0, 0.0, 0.5, 1.0, 0.0, Radius_geom),
                (0.5, 1.0, 0.0, 1.0, 1.0, 0.0, Radius_geom),
                (0.5, 0.0, 0.0, 0.5, 1.0, 0.0, Radius_geom),
                (0.0, 0.0, 1.0, 0.5, 0.0, 1.0, Radius_geom),
                (0.5, 0.0, 1.0, 1.0, 0.0, 1.0, Radius_geom),
                (0.0, 1.0, 1.0, 0.5, 1.0, 1.0, Radius_geom),
                (0.5, 1.0, 1.0, 1.0, 1.0, 1.0, Radius_geom),
                (0.5, 0.0, 1.0, 0.5, 1.0, 1.0, Radius_geom),
                (0.5, 0.0, 0.0, 0.5, 0.0, 1.0, Radius_geom),
                (0.5, 1.0, 0.0, 0.5, 1.0, 1.0, Radius_geom)]
    CubicV4 = [(0.5, 0.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
               (0.0, 0.5, 0.0, 0.5, 0.5, 0.0, Radius_geom),
               (0.5, 1.0, 0.0, 0.5, 0.5, 0.0, Radius_geom),
               (1.0, 0.5, 0.0, 0.5, 0.5, 0.0, Radius_geom),
               (0.5, 0.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
               (0.0, 0.5, 1.0, 0.5, 0.5, 1.0, Radius_geom),
               (0.5, 1.0, 1.0, 0.5, 0.5, 1.0, Radius_geom),
               (1.0, 0.5, 1.0, 0.5, 0.5, 1.0, Radius_geom),
               (0.5, 0.5, 0.0, 0.5, 0.5, 1.0, Radius_geom)]
    Newlattice = [
        (0.0, 0.0, 0.0, 0.25, 0.25, 0.25, Radius_geom),
        (0.25, 0.25, 0.25, 0.5, 0.0, 0.0, Radius_geom),
        (0.25, 0.25, 0.25, 0.0, 0.0, 0.5, Radius_geom),
        (0.25, 0.25, 0.25, 0.0, 0.5, 0.0, Radius_geom),
        (1.0, 0.0, 0.0, 0.75, 0.25, 0.25, Radius_geom),
        (0.75, 0.25, 0.25, 0.5, 0.0, 0.0, Radius_geom),
        (0.75, 0.25, 0.25, 1.0, 0.5, 0.0, Radius_geom),
        (0.75, 0.25, 0.25, 1.0, 0.0, 0.5, Radius_geom),
        (0.0, 1.0, 0.0, 0.25, 0.75, 0.25, Radius_geom),
        (0.25, 0.75, 0.25, 0.0, 0.5, 0.0, Radius_geom),
        (0.25, 0.75, 0.25, 0.5, 1.0, 0.0, Radius_geom),
        (0.25, 0.75, 0.25, 0.0, 1.0, 0.5, Radius_geom),
        (1.0, 1.0, 0.0, 0.75, 0.75, 0.25, Radius_geom),
        (0.75, 0.75, 0.25, 1.0, 0.5, 0.0, Radius_geom),
        (0.75, 0.75, 0.25, 0.5, 1.0, 0.0, Radius_geom),
        (0.75, 0.75, 0.25, 1.0, 1.0, 0.5, Radius_geom),
        (0.0, 0.0, 1.0, 0.25, 0.25, 0.75, Radius_geom),
        (0.25, 0.25, 0.75, 0.0, 0.5, 1.0, Radius_geom),
        (0.25, 0.25, 0.75, 0.5, 0.0, 1.0, Radius_geom),
        (0.25, 0.25, 0.75, 0.0, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 1.0, 0.75, 0.25, 0.75, Radius_geom),
        (0.75, 0.25, 0.75, 1.0, 0.0, 0.5, Radius_geom),
        (0.75, 0.25, 0.75, 0.5, 0.0, 1.0, Radius_geom),
        (0.75, 0.25, 0.75, 1.0, 0.5, 1.0, Radius_geom),
        (0.0, 1.0, 1.0, 0.25, 0.75, 0.75, Radius_geom),
        (0.25, 0.75, 0.75, 0.0, 1.0, 0.5, Radius_geom),
        (0.25, 0.75, 0.75, 0.5, 1.0, 1.0, Radius_geom),
        (0.25, 0.75, 0.75, 0.0, 0.5, 1.0, Radius_geom),
        (1.0, 1.0, 1.0, 0.75, 0.75, 0.75, Radius_geom),
        (0.75, 0.75, 0.75, 1.0, 1.0, 0.5, Radius_geom),
        (0.75, 0.75, 0.75, 0.5, 1.0, 1.0, Radius_geom),
        (0.75, 0.75, 0.75, 1.0, 0.5, 1.0, Radius_geom),
    ]
    Diamond = [
        (0.0, 0.0, 0.0, 0.25, 0.25, 0.25, Radius_geom),
        (0.25, 0.25, 0.25, 0.5, 0.5, 0.0, Radius_geom),
        (0.25, 0.25, 0.25, 0.0, 0.5, 0.5, Radius_geom),
        (0.25, 0.25, 0.25, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 0.0, 0.75, 0.25, 0.25, Radius_geom),
        (0.75, 0.25, 0.25, 0.5, 0.5, 0.0, Radius_geom),
        (0.75, 0.25, 0.25, 1.0, 0.5, 0.5, Radius_geom),
        (0.75, 0.25, 0.25, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 1.0, 0.0, 0.75, 0.75, 0.25, Radius_geom),
        (0.75, 0.75, 0.25, 0.5, 0.5, 0.0, Radius_geom),
        (0.75, 0.75, 0.25, 1.0, 0.5, 0.5, Radius_geom),
        (0.75, 0.75, 0.25, 0.5, 1.0, 0.5, Radius_geom),
        (0.0, 1.0, 0.0, 0.25, 0.75, 0.25, Radius_geom),
        (0.25, 0.75, 0.25, 0.5, 0.5, 0.0, Radius_geom),
        (0.25, 0.75, 0.25, 0.0, 0.5, 0.5, Radius_geom),
        (0.25, 0.75, 0.25, 0.5, 1.0, 0.5, Radius_geom),
        (0.0, 0.0, 1.0, 0.25, 0.25, 0.75, Radius_geom),
        (0.25, 0.25, 0.75, 0.5, 0.5, 1.0, Radius_geom),
        (0.25, 0.25, 0.75, 0.0, 0.5, 0.5, Radius_geom),
        (0.25, 0.25, 0.75, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 0.0, 1.0, 0.75, 0.25, 0.75, Radius_geom),
        (0.75, 0.25, 0.75, 0.5, 0.5, 1.0, Radius_geom),
        (0.75, 0.25, 0.75, 1.0, 0.5, 0.5, Radius_geom),
        (0.75, 0.25, 0.75, 0.5, 0.0, 0.5, Radius_geom),
        (1.0, 1.0, 1.0, 0.75, 0.75, 0.75, Radius_geom),
        (0.75, 0.75, 0.75, 0.5, 0.5, 1.0, Radius_geom),
        (0.75, 0.75, 0.75, 1.0, 0.5, 0.5, Radius_geom),
        (0.75, 0.75, 0.75, 0.5, 1.0, 0.5, Radius_geom),
        (0.0, 1.0, 1.0, 0.25, 0.75, 0.75, Radius_geom),
        (0.25, 0.75, 0.75, 0.5, 0.5, 1.0, Radius_geom),
        (0.25, 0.75, 0.75, 0.0, 0.5, 0.5, Radius_geom),
        (0.25, 0.75, 0.75, 0.5, 1.0, 0.5, Radius_geom)
    ]
    if (Lattice == 0):
        return BCC
    if (Lattice == 1):
        return Octet
    if (Lattice == 2):
        return OctetExt
    if (Lattice == 3):
        return OctetInt
    if (Lattice == 4):
        return BCCZ
    if (Lattice == 5):
        return Cubic
    if (Lattice == 6):
        return OctahedronZ
    if (Lattice == 7):
        return OctahedronZcross
    if (Lattice == 8):
        return Kelvin
    if (Lattice == 9):
        return CubicV2
    if (Lattice == 10):
        return CubicV3
    if (Lattice == 11):
        return CubicV4
    if (Lattice == 12):
        return Newlattice
    if (Lattice == 13):
        return Diamond

def Lattice_geometry_corrected(Lattice, Radius_geom):
    Lattice_margin = 0.0001
    LatGeom = Lattice_geometry(Lattice, Radius_geom)
    LatGeom_modifie = []
    for ligne in LatGeom:
        nouvelle_ligne = []
        for coord in ligne:
            if coord == 0.0:
                nouvelle_ligne.append(coord - Lattice_margin)
            elif coord == 1.0:
                nouvelle_ligne.append(coord + Lattice_margin)
            else:
                nouvelle_ligne.append(coord)
        LatGeom_modifie.append(tuple(nouvelle_ligne))
    return LatGeom_modifie

def Type_lattice(Lattice):
    if (Lattice == 0):
        Type = 'BCC'
    if (Lattice == 1):
        Type = 'Octet'
    if (Lattice == 2):
        Type = 'OctetExt'
    if (Lattice == 3):
        Type = 'OctetInt'
    if (Lattice == 4):
        Type = 'BCCZ'
    if (Lattice == 5):
        Type = 'Cubic'
    if (Lattice == 6):
        Type = 'OctahedronZ'
    if (Lattice == 7):
        Type = 'OctahedronZcross'
    if (Lattice == 8):
        Type = 'Kelvin'
    if (Lattice == 9):
        Type = 'CubicV2'
    if (Lattice == 10):
        Type = 'CubicV3'
    if (Lattice == 11):
        Type = 'CubicV4'
    if (Lattice == 12):
        Type = 'Newlattice'
    if (Lattice == 13):
        Type = 'Diamond'
    return Type

def getVectorOrientation(Lattice):
    if (Lattice == 0):
        VectorOrientation = [0,1,-1]
    elif (Lattice == 1):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 2):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 3):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 4):
        VectorOrientation = [0,1,-1]
    elif (Lattice == 5):
        VectorOrientation = [0,1,-1]
    elif (Lattice == 6):
        VectorOrientation = [1,-1,1]
    elif (Lattice == 7):
        VectorOrientation = [1,-1,-1]
    elif (Lattice == 8):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 9):
        VectorOrientation = [1,0,-1]
    elif (Lattice == 10):
        VectorOrientation = [1,0,-1]
    elif (Lattice == 11):
        VectorOrientation = [1,0,-1]
    else:
        VectorOrientation = [0,0,-1]
    return VectorOrientation

def GetCorrectionExteriorBeam(Lattice,TypeCorrection):
    if (Lattice == 1) or (Lattice == 2) or (Lattice == 5) or (Lattice == 7) or (Lattice == 8) or (Lattice == 10) or (Lattice == 11):
        if TypeCorrection == 1:
            return 1
        elif TypeCorrection == 2:
            return 2
        else:
            return 0
    else:
        return 0

def Color_lattice(Lattice):
    if (Lattice == 0):
        color = "#e60049"
    if (Lattice == 1):
        color = "#0bb4ff"
    if (Lattice == 2):
        color = "#50e991"
    if (Lattice == 3):
        color = "#ffa300"
    if (Lattice == 4):
        color = "#9b19f5"
    if (Lattice == 5):
        color = "#e6d800"
    if (Lattice == 6):
        color = "#32CD32"
    if (Lattice == 7):
        color = 'olivedrab'
    if (Lattice == 8):
        color = 'cadetblue'
    if (Lattice == 9):
        color = 'hotpink'
    if (Lattice == 10):
        color = 'lightcoral'
    if (Lattice == 11):
        color = 'darkslategray'
    if (Lattice == 12):
        color = 'blue'
    return color