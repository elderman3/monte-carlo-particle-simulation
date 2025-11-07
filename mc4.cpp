// surface and delta-tracking algorithm for transport
// Statistical estimate of neutron flux: track-length and collision
// External source simulation algorithm

/*
Demonstration:

    Cylinder
    - Hollow cylinder, air + water
    - Point neutron-source in center
    - Total leak rate + analog + implicit estimates of collisions + absorptions

    - Collision visualization
    - Criticality source simulation
    - Average neutron density



*/


// Surface tracking = effectively use neutron direction to find next collision surface and collide there.
// Delta tracking = take max S and use that to sample 
// These both need effectively 1 function inside simulation to do the interaction sampling. 

// Simulation is done either in batches (External source) or in cycles (Criticality source) [Removes neutrons if too many, dupes if few]

// Fix 2 old problems


// For stats:
/*
It is assumed that the simulated population size is divided into a number of equal size
batches. The statistical estimators (mean + standard deviation) are collected by averaging over the batch-wise results.
*/