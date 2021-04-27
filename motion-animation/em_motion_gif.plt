# Setting ----------------------------------------
reset
set term gif animate delay 5 size 1280,720
set output "em_motion.gif"
set margins 0, 0, 0, 0
unset key
unset grid

set xl "x" font"Times:Italic, 22"
set yl "y" font"Times:Italic, 22"

# Parameter ----------------------------------------
m = 2.0                         # [kg]
e = 0.1                         # [C]
k = 15                          # [N/Cm]
B = 65                          # [N/Am]
b = e*B/(2*m)
R1 = 10                         # Large radius [m]
R2 = 30                         # Small radius [m]
omg0 = sqrt(e*k/m)              # Angular velocity of z
omg1 = b + sqrt(b**2-omg0**2)
omg2 = b - sqrt(b**2-omg0**2)
r = 0.8                         # Radius of the ball
dt = 0.001                      # step
dh = dt/6.0
dec = 50                        # Decimation
mag = 0.15
limit = 30/dt                   # loop limit

DATA = "motion_data_gif.dat"
set print DATA

# Functions ----------------------------------------
# Lorentz force
f1(x, y, z, vx, vy, vz) = vx                # dx/dt
f2(x, y, z, vx, vy, vz) = vy                # dy/dt
f3(x, y, z, vx, vy, vz) = vz                # dz/dt
f4(x, y, z, vx, vy, vz) = e/m*(k*x-vy*B)    # dvx/dt
f5(x, y, z, vx, vy, vz) = e/m*(k*y+vx*B)    # dvy/dt
f6(x, y, z, vx, vy, vz) = -2*e*k*z          # dvz/dt

# Runge-Kutta 4th (Define rk_i(x, y, z, vx, vy, vz))
do for[i=1:6]{
    rki = "rk"
    fi  = "f".sprintf("%d", i)
    rki = rki.sprintf("%d(x, y, z, vx, vy, vz) = (\
        k1 = %s(x, y, z, vx, vy, vz),\
        k2 = %s(x + dt*k1/2., y + dt*k1/2., z + dt*k1/2., vx + dt*k1/2., vy + dt*k1/2., vz + dt*k1/2.),\
        k3 = %s(x + dt*k2/2., y + dt*k2/2., z + dt*k2/2., vx + dt*k2/2., vy + dt*k2/2., vz + dt*k2/2.),\
        k4 = %s(x + dt*k3, y + dt*k3, z + dt*k3, vx + dt*k3, vy + dt*k3, vz + dt*k3),\
        dh * (k1 + 2*k2 + 2*k3 + k4))", i, fi, fi, fi, fi)
    eval rki
}

# Time
Time(t) = sprintf("{/Times:Italic t} = %3.1f s", t)

# Plot ----------------------------------------
# Initial Value
t  = 0
x  = R1 + R2
y  = 0
z  = 20
vx = 0
vy = R1*omg1+R2*omg2
vz = 0

print x, y, z, vx, vy, vz            # Write initial value into DAT file

# Draw initiate state for 70 steps
do for [i = 1:70] {
    set label 1 Time(t) at screen 0.5, 0.90 font 'Times:Normal, 20'

    set multiplot
        # xyz-space
            set view 60, 30, 1, 1.2
            set origin 0.020,0.03
            set size 0.60, 0.60*1280/720
            set zl "z" font"Times:Italic, 22"
            set xyplane at -50

            # Lorentz force
            set arrow 1 from x, y, z to x+mag*e*k*x, y+mag*e*k*y, z-mag*2*e*k*z front lc rgb'red' lw 2
            set arrow 2 from x, y, z to x-mag*e*B*vy, y+mag*e*B*vx, z-mag*2*e*k*z front lc rgb 'blue' lw 2

            # Ball
            set obj 1 circ at x, y, z size r fc rgb 'black' fs solid front
            set obj 2 circ at x, y, -40 size r fc rgb 'black' fs solid front
            splot[-50:50][-50:50][-50:50] DATA using 1:2:3 with line linecolor rgb "black" lw 1

        # xy-plane
            set pm3d map
            set origin 0.55,0.11
            set size 0.48, 0.48*1280/720
            unset zl

            # Gray circles
            set obj 3 circ at 0, 0 size R2-R1 fc rgb 'gray50' fs transparent noborder behind
            set obj 4 circ at 0, 0 size R2    fc rgb 'gray50' fs transparent noborder behind
            set obj 5 circ at 0, 0 size R2+R1 fc rgb 'gray50' fs transparent noborder behind

            splot[-50:50][-50:50][-50:50] DATA using 1:2:3 with line linecolor rgb "black" lw 1
        unset multiplot

        do for[j=3:5]{
            unset obj j			# In order not to show obj 3,4,5 in xyz-space
        }
}

# Update for limit steps
do for [i = 1:limit] {
    t  = t  + dt
    x  = x  + rk1(x, y, z, vx, vy, vz)
    y  = y  + rk2(x, y, z, vx, vy, vz)
    z  = z  + rk3(x, y, z, vx, vy, vz)
    vx = vx + rk4(x, y, z, vx, vy, vz)
    vy = vy + rk5(x, y, z, vx, vy, vz)
    vz = vz + rk6(x, y, z, vx, vy, vz)

    print x, y, z, vx, vy, vz            # Write the values into DAT file

    if(i%dec==0){
        # Update time
        set label 1 Time(t) at screen 0.5, 0.93 font 'Times:Normal, 20'

        set multiplot
        # xyz-space
            set view 60, 30, 1, 1.2
            set origin 0.020,0.03
            set size 0.60, 0.60*1280/720
            set zl "z" font"Times:Italic, 22"
            set xyplane at -50

            # Lorentz force
            set arrow 1 from x, y, z to x+mag*e*k*x, y+mag*e*k*y, z-mag*2*e*k*z front lc rgb'red' lw 2
            set arrow 2 from x, y, z to x-mag*e*B*vy, y+mag*e*B*vx, z-mag*2*e*k*z front lc rgb 'blue' lw 2

            # Ball
            set obj 1 circ at x, y, z size r fc rgb 'black' fs solid front
            set obj 2 circ at x, y, -50 size r fc rgb 'black' fs solid front
            splot[-50:50][-50:50][-50:50] \
                 DATA using 1:2:3 with line linecolor rgb "black" lw 1, \
                 DATA using 1:2:(-50) with line linecolor rgb "black" lw 1

        # xy-plane
            set pm3d map
            set origin 0.55,0.11
            set size 0.48, 0.48*1280/720
            unset zl

            # Gray circles
            set obj 3 circ at 0, 0 size R2-R1 fc rgb 'gray50' fs transparent noborder behind
            set obj 4 circ at 0, 0 size R2    fc rgb 'gray50' fs transparent noborder behind
            set obj 5 circ at 0, 0 size R2+R1 fc rgb 'gray50' fs transparent noborder behind

            splot[-50:50][-50:50][-50:50] \
                 DATA using 1:2:3 with line linecolor rgb "black" lw 1
        unset multiplot

        do for[j=3:5]{
            unset obj j			# In order not to show obj 3,4,5 in xyz-space
        }
    }
}
set out