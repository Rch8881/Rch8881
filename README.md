### Satellite positioning demonstration 

# Satellite positioning
This is a demonstration on using MATLAB to decode the GPS ephemeris data, using the least-squares (LS) for the pseudorange measurements and finally resolve the user location.



## Introduction 
The raw signal received from satellite are put in two file namely 'eph.dat' and 'rcvr.dat'. This demonstration will try to resolve the user location by using these data. MATLAB will be used in this demonstration. 

## Instruction to use
1. Download the files

2. Make sure “eph.dat” and “rcvr.dat” are put in same folder named as ‘Data’ and place together with all others components including “assignment_v6.m”, “rot.m”,” wgslla2xyz.m”,” wgsxyz2enu.m” and” wgsxyz2lla.m”. 

3. Open the “assignment_v6.m” file in MATLAB

4. Press the "RUN" button

5. Results will be shown in the command window in MATLAB.

## Result
The locations for the user have been successfully resolved by the demonstration.

The starting position ECEF were -2694685.473m, -4293642.366m and 3857878.924m.

And the final solved position is -2700421.279m, -4292547.361m and 3855271.891m.

As the target position ECEF are -2700400.000m, -4292560.000m and 3855270.000m, the positioning error is 24.822m.









<!--
**Rch8881/Rch8881** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
