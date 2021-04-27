# Setting ----------------------------------------
reset
set term png size 1280, 720
set margins 0, 0, 0, 0
unset key
unset grid

set xl "x" font"Times:Italic, 22"
set yl "y" font"Times:Italic, 22"

# Parameter ----------------------------------------
m     = 2.0                     # [kg]
e     = 0.1                     # [C]
k     = 15                      # [N/Cm]
B     = 65                      # [N/Am]
b     = e*B/(2*m)
R1    = 10                      # Large [m]
R2    = 30                      # Small [m]
omg0  = sqrt(e*k/m)
omg1  = b + sqrt(b**2-omg0**2)
omg2  = b - sqrt(b**2-omg0**2)
mag = 0.15

array Name[2]  = ["E_field", "B_field"]
array Color[2] = ["red", "blue"]
array File[2]  = [".dat", ".png"]


# Functions ----------------------------------------
# Electric field
Ex(x, y, z) = -k*x
Ey(x, y, z) = -k*y
Ez(x, y, z) = 2*k*z

# Magnetic field
Bx(x, y, z) = 0
By(x, y, z) = 0
Bz(x, y, z) = B


# Plot ----------------------------------------
inc = 8        # increment
do for [f=1:2]{
    if(f==1){    # Electric field
        k = 15
        B = 0
        mag = 0.005
    } else {    # Magnetic field
        k = 0
        B = 65
        mag = 0.08
    }
    
    set print Name[f].File[1]
    
    do for [x=-40:40:inc] {
        do for[y=-40:40:inc] {
            do for[z=-40:40:inc]{                
                vx = mag*((f==1)?Ex(x, y, z):Bx(x, y, z))
                   vy = mag*((f==1)?Ey(x, y, z):By(x, y, z))
                vz = mag*((f==1)?Ez(x, y, z):Bz(x, y, z))
                print x, y, z, vx, vy, vz
            }
        }
    }
    
    set output Name[f].File[2]
    set multiplot
    # xyz-space
        set view 60, 30, 1, 1.2
        set origin 0.020,0.03
        set size 0.60, 0.60*1280/720
        set zl "z" font"Times:Italic, 22"
        set xyplane at -50
        
        splot[-50:50][-50:50][-50:50] \
             Name[f].File[1] using 1:2:3:4:5:6 w vec lc rgb Color[f] lw 1
        
    # xy-plane
        set pm3d map
        set origin 0.55,0.11
        set size 0.48, 0.48*1280/720
        unset zl
        
        splot[-50:50][-50:50][-50:50] \
             Name[f].File[1] using 1:2:3:4:5:6 w vec lc rgb Color[f] lw 1
    unset multiplot
}
set out    