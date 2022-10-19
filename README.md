# sensornetwork

This is an assignment from my Data Structures and Algorithms course in college. Below is the project description and objective. All code should be executable from any IDE and the main function call can be found in SensorNetwork.java and simply follow the prompts in terminal after running.

A window should appear that visualizes the network that you have generated after giving it the parameters needed in the terminal prompts. 

Description:

Network Model. The IoT sensor network is represented as a graph G(V,E), where V = {1, 2, ..., N} is the set of N nodes, and E is the set of edges. The sensor nodes are randomly generated in an x by y area. All the sensor nodes have the same transmission range Tr. That is, if the distance between two nodes are less than or equal to Tr, then they can communicate directly and are connected by an edge.

Out of N nodes p of them are randomly selected as storage-depleted data generating nodes, called data nodes, DN = {DN1, DN2, ..., DNp}, where DNi has si number of data packets. Each data is of 400 bytes (i.e., 3200 bits). We assume that each DN has zero storage capacity, while other sensor nodes, which are called storage nodes (SNs), have storage capacities – SNi has storage capacity mi, meaning it can store mi data packets. To preserve the data, all the data packets at the DNs must be offloaded to some SNs.

Energy Model. All the sensor nodes are powered by battery, which has limited amount of energy reservoir. Sensor nodes communicate with each other using wireless communication, which costs the battery power of sensor nodes. For two sensor nodes i an j who are within each other’s transmission range, they can communicate by sending and receiving data packets. Assume their distance is l meters. When node i sends k-bit data packet to its neighbor node j over their l-meter distance, both i and j spend their battery power as follows:

The transmission energy spent by node i is ET (k, l ) = Eelec  k +Eamp  k  l2 and the receiving energy spent by node j is ER (k ) = Eelec  k. Here, Eamp = 100 pJ/bit/m2 is the energy consumption of one bit on i’s transmit amplifier and Eelec = 100nJ/bit is the
energy consumption of one bit on i’s transmitter circuit and j’s receiver circuit.

Therefore for any edge (i, j) with length l in the sensor network, its weight (i.e., its cost) = ET (k, l ) + ER (k) = 2  Eelec  k +Eamp  k  l2.

Programming Objective: When executed, your program should prompt to ask
1. the width x and length y of the sensor network (e.x., 50meter x 50meter)
2. number of sensor nodes: N
3. transmission range in meters: Tr
4. number of DNs: p, and number of data packets each DN has: q. Here we assume each DN has the same number of data packets.
5. Storage capacity of each storage node: m. Here we assume each storage node has the same storage capacity.

First, your program still needs to check if the sensor network graph is connected (i.e.. connectivity). If not, it prints out a message “the network is not connected”, and asks the user to input again.
Second, your program should check if the user inputs satisfy: p  q <= (N-p)  m (i.e., feasibility). That is, there should have enough storage spaces in the entire network to store all the data packets generated. If not, it prints out a message “there is not enough storage in the network’, and asks the user to input again.
After both connectivity and feasibility are satisfied, next, your program should list the IDs of DNs and the IDs of storage nodes. It then asks user to input a DN node and a storage node as follows:

6. Please input the IDs of a DN and a SN:

7. Please choose algorithm to execute: 0 for Dijkstra’s shortest path algorithm, 1 for Bellman-Ford dynamic programming algorithm, and 2 and for finding a shortest path between them with k edges. If input is 2, then ask to input k.
Your program should output the minimum-energy data offloading path; energy cost of offloading one data packet from the DN to the SN along this path, and the total energy cost of offloading all this DN’s q data packets to the SN. Note the energy cost along a path is the sum of the weights of all the edges of the path.
