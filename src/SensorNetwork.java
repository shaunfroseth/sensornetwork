/**
 * SensorNetwork.java
 * Updated by Amado Rodriguez III, Shaun Froseth, Nelly Sanchez-Cruz
 * Ver. 3.0 for proj 2 4/15/21
 * Ver. 2.0 for proj 1 3/19/21
 * Ver. 1.0 for proj 0 3/4/21
 * CSC 401
 */

import java.awt.Dimension;
//import java.sql.SQLOutput;
import java.util.*;

public class SensorNetwork {
	
	private Map<Integer, Axis> nodes = new LinkedHashMap<Integer, Axis>();
	Set<Integer> dn = new HashSet<>();
	Map<Integer, Boolean> discovered = new HashMap<Integer, Boolean>();
	Map<Integer, Boolean> explored = new HashMap<Integer, Boolean>();
	Map<Integer, Integer> parent = new HashMap<Integer, Integer>();
	Map<Integer, Integer> path = new HashMap<>();
	Stack<Integer> s = new Stack<Integer>();
	PriorityQueue priorityQueue;
	List<List<Node>> nodeAdjList;
	Set<Integer> visited;
	double[] dist;
	static int traversal;

	public static void main(String[] args) {
		Scanner scan = new Scanner(System.in);
		boolean connected = false;
		boolean feasible = false;
		SensorNetwork sensor = new SensorNetwork();
		Map<Integer, Set<Integer>> adjacencyList1 = new LinkedHashMap<>();
		double width, height;
		int numberOfNodes = 0, sourceDN = 0, p = 0, q = 0;

		while(!connected || !feasible) {
			System.out.println("Enter the width (in meters):");
			width = scan.nextDouble();

			System.out.println("Enter the height (in meters):");
			height = scan.nextDouble();

			System.out.println("Enter the number of nodes:");
			numberOfNodes = scan.nextInt();

			System.out.println("Enter the Transmission range (in meters):");
			int transmissionRange = scan.nextInt();

			// New to Proj 1
			System.out.println("Enter the number of data nodes:");
			p = scan.nextInt();

			System.out.println("Enter the number of data packets each data node has:");
			q = scan.nextInt();

			System.out.println("Enter the storage capacity of each storage node:");
			int m = scan.nextInt();

			System.out.println("Enter the graph traversal technique:");
			System.out.println("Recursive DFS: 0");
			System.out.println("DFS Using Stack: 1");
			System.out.println("BFS Using Queue: 2");
			traversal = scan.nextInt();

			sensor.populateNodes(numberOfNodes, width, height);

			System.out.println("\nNode List:");
			for (int key : sensor.nodes.keySet()) {
				Axis ax = sensor.nodes.get(key);
				System.out.println("Node:" + key + ", xAxis: " + ax.getxAxis() + ", yAxis: " + ax.getyAxis());
			}

			sensor.populateAdjacencyList(numberOfNodes, transmissionRange, adjacencyList1);

			connected = sensor.executeAlgorithms(width, height, adjacencyList1);
			feasible = sensor.feasibility(numberOfNodes, p, q, m);
			if(!feasible)
				System.out.println("There is not enough storage in the network. Please try again.");
		}

		sensor.generateRandomDNs(numberOfNodes, p);
		System.out.println("Data Nodes: ");
		for (int d : sensor.dn) {
			System.out.println(d);
		}
		System.out.println("Storage Nodes: ");
		for (int s = 1; s <= numberOfNodes; s++) {
			if (!sensor.dn.contains(s)) {
				System.out.println(s);
			}
		}

		/*
		System.out.println("\nAdjacency List: ");
		for (int i : adjacencyList1.keySet()) {
			System.out.print(i);
			if (!adjacencyList1.isEmpty()) {
				for (int j : adjacencyList1.get(i)) {
					System.out.print("->" + j);
				}
			}
			System.out.println();
		}
		*/

		// array of distances
		sensor.dist = new double[numberOfNodes + 1];
		sensor.visited = new HashSet<>();
		sensor.priorityQueue = new PriorityQueue(numberOfNodes);

		// Copying adjacency list map to a new list of lists of nodes
		sensor.nodeAdjList = new ArrayList<>();
		sensor.nodeAdjList.add(null);
		for (Integer key : adjacencyList1.keySet()) {
			List<Node> nodeList = new ArrayList<>();
			nodeList.add(new Node(0, 0));
			for (Integer node : adjacencyList1.get(key)){
				nodeList.add(new Node(node, 0));
			}
			sensor.nodeAdjList.add(nodeList);
		}

		System.out.println("Please input the ID of a DN: ");
		sourceDN = scan.nextInt();
		System.out.println("Please input the ID of a target SN: ");
		int targetSN = scan.nextInt();

		System.out.println("Enter the shortest path algorithm desired:");
		System.out.println("Dijkstra's shortest path: 0");
		System.out.println("Bellman-Ford dynamic programming: 1");
		System.out.println("Shortest path between them with k edges: 2");
		int algo = scan.nextInt();

		// New to project 2

		switch (algo) {
			case 0:
				sensor.dijkstraAlgorithm(sourceDN, numberOfNodes);
				break;
			case 1:
				sensor.bellmanFord(sourceDN, numberOfNodes);
				break;
			case 2:
				System.out.println("Please enter the number of edges k: ");
				int k = scan.nextInt();

				double weight = sensor.kEdges(sourceDN, numberOfNodes, targetSN, k);
				double cost = 0.00064 + 0.00000032 * weight * weight;
				if(weight == Integer.MAX_VALUE){
					System.out.println("There is no path with " + k + " edges.");
				}
				else {
					System.out.print("The cost of offloading one data packet from DN# ");
					System.out.println(sourceDN + " to SN# " + targetSN + " using " + k + " edges is: ");
					System.out.println(cost + " Joules. ");
					System.out.println("The minimum energy cost of offloading all data packets from this node is:");
					cost *= q;
					System.out.println(cost + " Joules. ");
				}
		}

		scan.close();

		if(algo != 2) {
			double cost = 0.00064 + 0.00000032 * sensor.dist[targetSN] * sensor.dist[targetSN];
			System.out.print("The minimum energy cost of offloading one data packet from DN# ");
			System.out.println(sourceDN + " to SN# " + targetSN + " is: ");
			System.out.println(cost + " Joules. ");
			System.out.println("The minimum energy cost of offloading all data packets from this node is:");
			cost *= q;
			System.out.println(cost + " Joules. ");

			System.out.println("Using this path: ");
			System.out.print(targetSN + " -> ");
			int root = sensor.path.get(targetSN);
			while (root != sourceDN) {
				System.out.print(root + " -> ");
				root = sensor.path.get(root);
			}
			System.out.println(root);
		}
	}

	boolean executeAlgorithms(double width, double height, Map<Integer, Set<Integer>> adjList) {

		switch (traversal) {
			case 0:
				System.out.println("\nExecuting Recursive DFS Algorithm");
				break;
			case 1:
				System.out.println("\nExecuting DFS Algorithm with Stack Algorithm");
				break;
			case 2:
				System.out.println("\nExecuting BFS with Queue Algorithm");
				break;
			default:
				System.out.println("\nNo valid graph traversal technique selected");
				break;
		}
		List<Set<Integer>> connectedNodes = new ArrayList<Set<Integer>>();
		for(int node: adjList.keySet()) {
			Set<Integer> connectedNode = new LinkedHashSet<Integer>();

			// Switch to check connections based on algorithm selected
			switch (traversal) {
				case 0:
					recursiveDFS(node, connectedNode, adjList);
					break;
				case 1:
					stackDFS(node, connectedNode, adjList);
					break;
				case 2:
					queueBFS(node, connectedNode, adjList);
					break;
			}

			if(!connectedNode.isEmpty()) {
				connectedNodes.add(connectedNode);
			}
		}

		drawGraph(width, height, adjList);

		if(connectedNodes.size() == 1) {
			System.out.println("The network is fully connected with one connected component.");
		}
		else {
			System.out.println("The network is not connected. Please input again.");
			return false;
		}

		// System.out.println("There are " + connectedNodes.size() + " connected components");
		for(Set<Integer> list: connectedNodes) {
			System.out.println(list);
		}

		//Draw sensor network graph

		return true;
	}

	void recursiveDFS(int u, Set<Integer> connectedNode, Map<Integer, Set<Integer>> adjList) {
		
		if(!s.contains(u) && !explored.containsKey(u)) {
			s.add(u);
			discovered.put(u, true);
		} 
		
		while(!s.isEmpty()) {
			if(!explored.containsKey(u)) {
				List<Integer> list = new ArrayList<Integer>(adjList.get(u));
				for(int v: list) {

					if(!discovered.containsKey(v)) {
						s.add(v);
						discovered.put(v, true);

						parent.putIfAbsent(v, u);
						recursiveDFS(v, connectedNode, adjList);
					}
					else if(list.get(list.size()-1) == v) {
						if( parent.containsKey(u)) {
							explored.put(u, true);
							s.removeElement(u);

							connectedNode.add(u);
							recursiveDFS(parent.get(u), connectedNode, adjList);
						}
					}
				}
				if(!explored.containsKey(u))
					explored.put(u, true);
				s.removeElement(u);
				connectedNode.add(u);
			}
		}
			
	}

	/**
	 * stackDFS implements a stack to iteratively perform a Depth First Search using the current node as the root.
	 * @param u the current node
	 * @param connectedNode a set of connected nodes updated with the current node
	 * @param adjList an adjacency list of all nodes in the graph
	 */
	void stackDFS(int u, Set<Integer> connectedNode, Map<Integer, Set<Integer>> adjList){
		if(explored.containsKey(u)){
			return;
		}

		s.push(u);
		int currentNode = u;

		while(!s.isEmpty()){
			List<Integer> list = new ArrayList<>(adjList.get(currentNode));
			for(int v: list){
				if(!explored.containsKey(v)){
					s.push(v);
				}
			}
			connectedNode.add(currentNode);
			explored.put(currentNode, true);
			currentNode = s.pop();
		}
	}

	/**
	 * queueBFS implements a queue to iteratively perform a Breadth First Search using the current node as the root.
	 * @param u the current node
	 * @param connectedNode a set of connected nodes updated with the current node
	 * @param adjList an adjacency list of all nodes in the graph
	 */
	void queueBFS(int u, Set<Integer> connectedNode, Map<Integer, Set<Integer>> adjList){
		if(explored.containsKey(u)){
			return;
		}

		Queue<Integer> q = new LinkedList<>();
		q.add(u);
		int currentNode = u;

		while(!q.isEmpty()){
			List<Integer> list = new ArrayList<>(adjList.get(currentNode));
			for(int v: list){
				if(!explored.containsKey(v)){
					q.add(v);
				}
			}
			connectedNode.add(currentNode);
			explored.put(currentNode, true);
			currentNode = q.poll();
		}
	}

	void populateNodes(int nodeCount, double width, double height) {
		Random random = new Random();
		
		for(int i = 1; i <= nodeCount; i++) {
			Axis axis = new Axis();
			int scale = (int) Math.pow(10, 1);
			double xAxis =(0 + random.nextDouble() * (width - 0));
			double yAxis = 0 + random.nextDouble() * (height - 0);
			
			xAxis = (double)Math.floor(xAxis * scale) / scale;
			yAxis = (double)Math.floor(yAxis * scale) / scale;
			
			axis.setxAxis(xAxis);
			axis.setyAxis(yAxis);
			
			nodes.put(i, axis);	
		}
	}

	/**
	 * Updates the set dn with p random unique Integers from the pool of n nodes.
	 * Each Integer in dn represents the ID of an individual node.
	 * @param n the total number of nodes in the network
	 * @param p the number of Data Nodes in the network
	 */
	void generateRandomDNs(int n, int p) {
		List<Integer> list = new ArrayList<>();
		for (int i = 1; i <= n; i++) {
			list.add(i);
		}
		Collections.shuffle(list);
		for (int i = 0; i < p; i++) {
			dn.add(list.get(i));
		}
	}

	/**
	 * Checks to see if the current values of n, p, q, and m are feasible; ie. there are
	 * enough storage spaces in the network to store all the data packets generated.
	 * @param n the total number of nodes in the network
	 * @param p the total number of Data Nodes (DN) in the network
	 * @param q the number of data packets each DN has
	 * @param m the storage capacity of each storage node
	 * @return true iff feasible
	 */
	boolean feasibility(int n, int p, int q, int m) {
		return p * q <= (n - p) * m;
	}
	
	void populateAdjacencyList(int nodeCount, int tr, Map<Integer, Set<Integer>> adjList) {
		for(int i=1; i<= nodeCount; i++) {
			adjList.put(i, new HashSet<>());
		}
		
		for(int node1: nodes.keySet()) {
			Axis axis1 = nodes.get(node1);
			for(int node2: nodes.keySet()) {
				Axis axis2 = nodes.get(node2);
				
				if(node1 == node2) {
					continue;
				}
				double xAxis1 = axis1.getxAxis();
				double yAxis1 = axis1.getyAxis();
					
				double xAxis2 = axis2.getxAxis();
				double yAxis2 = axis2.getyAxis();
				
				double distance =  Math.sqrt(((xAxis1-xAxis2)*(xAxis1-xAxis2)) + ((yAxis1-yAxis2)*(yAxis1-yAxis2)));
				
				if(distance <= tr) {
					Set<Integer> tempList = adjList.get(node1);
					tempList.add(node2);
					adjList.put(node1, tempList);
						
					tempList = adjList.get(node2);
					tempList.add(node1);
					adjList.put(node2, tempList);
				}
			}
		}
	}

	/**
	 * drawGraph Draws a graph! :)
	 * @param width width of the graph
	 * @param height height of the graph
	 * @param adjList an adjacency list of all nodes in the graph
	 */
	void drawGraph(double width, double height, Map<Integer, Set<Integer>> adjList) {
		SensorNetworkGraph graph = new SensorNetworkGraph();
		graph.setGraphWidth(width);
		graph.setGraphHeight(height);
		graph.setNodes(nodes);
		graph.setAdjList(adjList);
		graph.setPreferredSize(new Dimension(960, 800));
		Thread graphThread = new Thread(graph);
		graphThread.start();
	}

	/**
	 * Dijkstra Algorithm for finding the shortest path between nodes in the network
	 * @param dn The source data node whose shortest path to the other nodes you are interested in
	 * @param n the number of nodes in the network
	 */
	void dijkstraAlgorithm(int dn, int n) {
		for (int i = 0; i <= n; i++) {
			dist[i] = Integer.MAX_VALUE;
		}

		priorityQueue.insert(new Node(dn, 0));
		dist[dn] = 0;

		while(visited.size() != n) {
			int root = priorityQueue.remove().node;
			visited.add(root);
			dijkstraAdjNodes(root);
		}
	}

	/**
	 * A helper method used by dijkstraAlgorithm that updates dist[] with shortest distances and the priority queue
	 * with new Nodes
	 * @param root The current Node to be processed
	 */
	void dijkstraAdjNodes(int root) {
		for (int i = 1; i < nodeAdjList.get(root).size(); i++) {
			Node vertex = nodeAdjList.get(root).get(i);

			if(!visited.contains(vertex.node)) {
				Axis axis1 = nodes.get(root);
				Axis axis2 = nodes.get(vertex.node);

				double xAxis1 = axis1.getxAxis();
				double yAxis1 = axis1.getyAxis();

				double xAxis2 = axis2.getxAxis();
				double yAxis2 = axis2.getyAxis();

				double distance =  Math.sqrt(((xAxis1-xAxis2)*(xAxis1-xAxis2)) + ((yAxis1-yAxis2)*(yAxis1-yAxis2)));

				double newdist = dist[root] + distance;

				if(newdist < dist[vertex.node]) {
					dist[vertex.node] = newdist;
					path.put(vertex.node, root);
				}

				priorityQueue.insert(new Node(vertex.node, dist[vertex.node]));
				visited.add(vertex.node);
			}
		}
	}

	/**
	 * Uses the Bellman-Ford algorithm to find the shortest distance from dn to all other vertices.
	 * @param dn The source data node.
	 * @param n Number of nodes in the graph.
	 */
	void bellmanFord(int dn, int n){
		for (int i = 0; i <= n; i++) {
			dist[i] = Integer.MAX_VALUE;
		}
		// Distance to itself should be 0
		dist[dn] = 0;

		// Loop through each Node's adjacency list
		for(int i = 1; i < nodeAdjList.size(); i++) {
			for(int j = 1; j < nodeAdjList.get(i).size(); j++) {
				Node vertex = nodeAdjList.get(i).get(j);
				if(vertex == null)
					continue;

				Axis axis1 = nodes.get(i);
				Axis axis2 = nodes.get(vertex.node);

				double xAxis1 = axis1.getxAxis();
				double yAxis1 = axis1.getyAxis();

				double xAxis2 = axis2.getxAxis();
				double yAxis2 = axis2.getyAxis();

				double distance =  Math.sqrt(((xAxis1-xAxis2)*(xAxis1-xAxis2)) + ((yAxis1-yAxis2)*(yAxis1-yAxis2)));

				if (dist[i] != Integer.MAX_VALUE && dist[i] + distance < dist[vertex.node]) {
					dist[vertex.node] = dist[i] + distance;
					path.put(vertex.node, i);
				}
			}
		}
	}

	/**
	 * Uses a dynamic programming algorithm to find the shortest distance from dn to the target node t.
	 * @param dn The source data node.
	 * @param n Number of nodes in the graph.
	 * @param t The target node.
	 * @param k The number of edges required.
	 */
	double kEdges(int dn, int n, int t, int k) {
		double[][] dynamicDistance = new double[n+1][k+1];
		dynamicDistance[t][0] = 0;
		for (int i = 0; i < n + 1; i++) {
			if(i != t)
				dynamicDistance[i][0] = Integer.MAX_VALUE;
		}
		for (int i = 1; i < nodeAdjList.get(t).size(); i++) {
			if(t == i) {
				continue;
			}
			Node vertex = nodeAdjList.get(t).get(i);
			if(vertex != null) {
				Axis axis1 = nodes.get(t);
				Axis axis2 = nodes.get(vertex.node);

				double xAxis1 = axis1.getxAxis();
				double yAxis1 = axis1.getyAxis();

				double xAxis2 = axis2.getxAxis();
				double yAxis2 = axis2.getyAxis();

				double distance = Math.sqrt(((xAxis1 - xAxis2) * (xAxis1 - xAxis2)) + ((yAxis1 - yAxis2) * (yAxis1 - yAxis2)));
				dynamicDistance[i][1] = distance;
			}
			else {
				dynamicDistance[i][1] = Integer.MAX_VALUE;
			}
		}

		for (int edges = 1; edges < k + 1; edges++) {
			for (int currentNode = 1; currentNode < n + 1; currentNode++) {
				dynamicDistance[currentNode][edges] = Integer.MAX_VALUE;
				for (Node adjNode: nodeAdjList.get(currentNode)) {
					if(adjNode == null){
						continue;
					}
					if(dynamicDistance[adjNode.node][edges - 1] != Integer.MAX_VALUE) {
						double weight = dynamicDistance[adjNode.node][edges - 1];

						Axis axis1 = nodes.get(currentNode);
						Axis axis2 = nodes.get(adjNode.node);

						double xAxis1 = axis1.getxAxis();
						double yAxis1 = axis1.getyAxis();

						double xAxis2 = axis2.getxAxis();
						double yAxis2 = axis2.getyAxis();

						double distance =  Math.sqrt(((xAxis1-xAxis2)*(xAxis1-xAxis2)) + ((yAxis1-yAxis2)*(yAxis1-yAxis2)));

						dynamicDistance[currentNode][edges] = Math.min(dynamicDistance[currentNode][edges], weight + distance);
					}
				}
			}
		}
		return dynamicDistance[dn][k];
	}

}
