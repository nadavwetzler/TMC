# The Time Magnitude Clustering (TMC) method 
The method gathers the clusters using a linkage algorithm prescribed by magnitude-dependent spatial and temporal distances following the popular Gardner and Knopoff (1974) declustering method. 
Rather than using the original specific spatial and temporal conditions, tailored to Southern California, we retain the basic approach of the Gardner and Knopoff formulation, while modifying the conditions for spatial and temporal distances. 
The spatial distance for triggering between earthquake pairs is based on 
a radius from the mainshock hypocenter given by a scaled value of the magnitude-based mainshock rupture length (LWnC) estimated by Wells and Coppersmith (1994), adjusted to account for the huge rupture lengths of some megathrust earthquakes:

R=LWnC*q           q=1.0, if M≥9; 1.5, if 8≥M>9; 2.0, if M<8 


, and the temporal triggering window is defined in Equation 2. In that sense, clusters are defined in a cumulative fashion, where the final duration and geographic size of each cluster is independent of a fixed spatial or time window, and do not associate exclusively with a specific mainshock. The temporal window length is defined as: 

t=exp(M)  

where t is the time window (in days). This results in triggering windows of 148, 403, 1096, 2980, and 8103 days for the magnitude range 5 – 9, 
all exceeding the 60-day window for the fixed-time procedure. 

![TMC](https://user-images.githubusercontent.com/88764899/162423850-9ec5007c-4561-45f5-a335-c99db70098f3.png)
