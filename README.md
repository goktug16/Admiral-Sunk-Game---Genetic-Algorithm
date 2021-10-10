# Admiral-Sunk-Game---Genetic-Algorithm


In the basic logic of the admiral sunk game, first of all, Structures such as enemy and ally ships, as well as bridges were created. Our goal is to find the best firing order to take down enemy ships and bridges. But we do not want to sink ally ships so there is a point system where we got points for taking down enemy ships and bridges. Bridges count as strategic places for supply transportation.

The 9*9 matrix given to us is thought of as a sea. In this sea, there may be some type of vehicles or bridges for each coordinate, these are in order;
1 –  Ally ships
2 – Enemy ships
3 – Bridges

These targets in the sea are wanted to be sunk by the artillery unit. Artillery fire starts from a certain point and from that point, it is aimed to hit the enemy ships and bridges. 
The gunner aims to hit enemy ships and bridges while avoiding hitting friendly ships. While the cannonballs continue, they cannot make a sharp turn to another point, as they are started from a certain point.

This situation has been created considering that it will be difficult to rotate the cannon. It is crucial to gain time between each shot. The angle of the cannon needs to change at a minimum after each shot.
Considering these circumstances, a possible cannon shooting pattern is found using genetic algorithm and this shooting pattern is applied.
