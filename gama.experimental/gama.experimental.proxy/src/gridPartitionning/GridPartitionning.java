package gridPartitionning;

import java.util.ArrayList;
import java.util.List;

import gama.annotations.precompiler.IConcept;
import gama.annotations.precompiler.IOperatorCategory;
import gama.core.runtime.IScope;
import gama.annotations.precompiler.GamlAnnotations.doc;
import gama.annotations.precompiler.GamlAnnotations.operator;
import gama.gaml.types.IType;
import java.util.*;

public class GridPartitionning {

	static Random random;
	
	class ClusterStats {
	    final int minClusterSize;
	    final int maxClusterSize;
	    final double avgClusterSize;
	    final double stdDevClusterSize;
	    final int borderCellCount;
	    final Map<Integer, Integer> borderCellsPerCluster;
	    final Map<Integer, Integer> sizesPerCluster;
	    final int totalCells;
	    final double borderPct;
	    final int numClusters;
	    final double sizeVariationCoeff;  // Coefficient of variation (stdDev/mean)

	    ClusterStats(int min, int max, double avg, double stdDev, int borders, int adj, 
	                 Map<Integer, Integer> borderPerCluster,
	                 Map<Integer, Integer> sizesPerCluster,
	                 int totalCells, double borderPct) {
	        this.minClusterSize = min;
	        this.maxClusterSize = max;
	        this.avgClusterSize = avg;
	        this.stdDevClusterSize = stdDev;
	        this.borderCellCount = adj;
	        this.borderCellsPerCluster = borderPerCluster;
	        this.sizesPerCluster = sizesPerCluster;
	        this.totalCells = totalCells;
	        this.borderPct = borderPct;
	        this.numClusters = sizesPerCluster.size();
	        this.sizeVariationCoeff = avg != 0 ? (stdDev/avg) : 0;
	    }

	    @Override
	    public String toString() {
	        StringBuilder sb = new StringBuilder();
	        sb.append(String.format(
	            "Cluster Stats:\n" +
	            "Total cells: %d\n" +
	            "Number of clusters: %d\n" +
	            "Min cells per cluster: %d\n" +
	            "Max cells per cluster: %d\n" +
	            "Avg cells per cluster: %.2f\n" +
	            "StdDev of cluster sizes: %.2f\n" +
	            "Size variation coefficient: %.2f\n" +
	            "Border Cells: %d (%.2f%%)\n\n",
	            totalCells, numClusters, minClusterSize, maxClusterSize, 
	            avgClusterSize, stdDevClusterSize, sizeVariationCoeff,
	            borderCellCount, borderPct));

	        // Add detailed size distribution
	        sb.append("Size distribution:\n");
	        sizesPerCluster.entrySet().stream()
	            .sorted(Map.Entry.comparingByKey())
	            .forEach(e -> sb.append(String.format("Cluster %d: %d cells (%.1f%% of total)\n", 
	                e.getKey(), e.getValue(), (100.0 * e.getValue() / totalCells))));
	        
	        return sb.toString();
	    }
	}

	static public class GRID {
	    private Integer[][] grid; 
	    private int rows;
	    private int cols;
	    private Integer[][] clusterGrid;
	    private Integer[][] scoreGrid; // Add score grid for CAMMISOL
	    private List<Point> borderCells = new ArrayList<>();
	    private List<Point> adjCell = new ArrayList<>();
	    private Map<Point, Set<Integer>> borderCellClusters = new HashMap<>();
	    

	    public GRID(int rows, int cols) {
	        this.rows = rows;
	        this.cols = cols;
	        this.grid = new Integer[rows][cols];
	        this.scoreGrid = new Integer[rows][cols]; // Initialize score grid
	        initializeGrid();
	    }
	    
	    // Inner class for agent cells with scores
	    private static class AgentCell {
	        int row, col, score;
	        AgentCell(int row, int col, int score) {
	            this.row = row;
	            this.col = col;
	            this.score = score;
	        }
	    }
		
	    private static class Point {
	        int row, col;
	        Point(int row, int col) {
	            this.row = row;
	            this.col = col;
	        }
	        double distanceTo(Point other) {
	            return Math.sqrt(Math.pow(row - other.row, 2) + Math.pow(col - other.col, 2));
	        }
	    }

	    private void initializeGrid() {
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                grid[i][j] = 0;
	                scoreGrid[i][j] = 0; // Initialize scores to 0
	            }
	        }
	    }

	    private void initializeClusters(int K) {
	        // Initialize clusterGrid with -1
	        clusterGrid = new Integer[rows][cols];
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                clusterGrid[i][j] = -1;
	            }
	        }
	    }

	    public Map<Integer, Integer> getClusterSizes() {
	        Map<Integer, Integer> sizes = new HashMap<>();
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (grid[i][j] == 0 && clusterGrid[i][j] != -1) {
	                    sizes.merge(clusterGrid[i][j], 1, Integer::sum);
	                }
	            }
	        }
	        return sizes;
	    }

	    int getTotalCells()
	    {
	        return rows * cols;
	    }
	    
	    // Methods for score management
	    public void setScore(int row, int col, int score) {
	        if (row >= 0 && row < rows && col >= 0 && col < cols) {
	            scoreGrid[row][col] = score;
	        }
	    }
	    
	    public int getScore(int row, int col) {
	        if (row >= 0 && row < rows && col >= 0 && col < cols) {
	            return scoreGrid[row][col];
	        }
	        return 0;
	    }
	    
	    public void setScoreGrid(Integer[][] scores) {
	        if (scores.length == rows && scores[0].length == cols) {
	            for (int i = 0; i < rows; i++) {
	                for (int j = 0; j < cols; j++) {
	                    scoreGrid[i][j] = scores[i][j];
	                }
	            }
	        }
	    }
	    
	    public Integer[][] getScoreGrid() {
	        Integer[][] copy = new Integer[rows][cols];
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                copy[i][j] = scoreGrid[i][j];
	            }
	        }
	        return copy;
	    }

	    public void divideGridToKMEANClusters(int K) {
	        List<Point> allCells = new ArrayList<>();
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                allCells.add(new Point(i, j));
	            }
	        }

	        if (allCells.isEmpty()) return;

	        // Initialize cluster grid
	        initializeClusters(K);

	        // Initialize K random centroids from all cells
	        List<Point> centroids = new ArrayList<>();
	        Set<Integer> usedIndices = new HashSet<>();
	        for (int i = 0; i < K && i < allCells.size(); i++) {
	            int index;
	            do {
	                index = random.nextInt(allCells.size());
	            } while (!usedIndices.add(index));
	            centroids.add(allCells.get(index));
	        }

	        // K-means iteration
	        boolean changed;
	        int maxIterations = 100;
	        do {
	            changed = false;
	            // Assign each cell to nearest centroid
	            Map<Integer, List<Point>> clusterPoints = new HashMap<>();
	            for (Point cell : allCells) {
	                int nearestCluster = 0;
	                double minDistance = Double.MAX_VALUE;
	                
	                for (int i = 0; i < centroids.size(); i++) {
	                    double distance = cell.distanceTo(centroids.get(i));
	                    if (distance < minDistance) {
	                        minDistance = distance;
	                        nearestCluster = i;
	                    }
	                }

	                clusterPoints.computeIfAbsent(nearestCluster, k -> new ArrayList<>()).add(cell);
	                int oldCluster = clusterGrid[cell.row][cell.col];
	                if (oldCluster != nearestCluster) {
	                    clusterGrid[cell.row][cell.col] = nearestCluster;
	                    changed = true;
	                }
	            }

	            // Update centroids
	            for (int i = 0; i < K; i++) {
	                List<Point> clusterCells = clusterPoints.getOrDefault(i, new ArrayList<>());
	                if (!clusterCells.isEmpty()) {
	                    double avgRow = clusterCells.stream().mapToInt(p -> p.row).average().orElse(0);
	                    double avgCol = clusterCells.stream().mapToInt(p -> p.col).average().orElse(0);
	                    centroids.set(i, new Point((int)avgRow, (int)avgCol));
	                }
	            }

	            maxIterations--;
	        } while (changed && maxIterations > 0);

	        // Count active clusters and remap indices
	        Set<Integer> activeClusters = new HashSet<>();
	        for (Point p : allCells) {
	            if (clusterGrid[p.row][p.col] != -1) {
	                activeClusters.add(clusterGrid[p.row][p.col]);
	            }
	        }
	        
	        // Create index mapping
	        Map<Integer, Integer> indexMap = new HashMap<>();
	        int newIndex = 0;
	        for (int oldIndex : activeClusters) {
	            indexMap.put(oldIndex, newIndex++);
	        }

	        // Remap cluster indices
	        for (Point p : allCells) {
	            if (clusterGrid[p.row][p.col] != -1) {
	                clusterGrid[p.row][p.col] = indexMap.get(clusterGrid[p.row][p.col]);
	            }
	        }

	        // After remapping indices, calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (grid[i][j] == 0 && clusterGrid[i][j] != -1) {  // Only consider cells in clusters
	                    if (isBorderCell(i, j, false)) {
	                        Point p = new Point(i, j);
	                        adjCell.add(p);
	                    }
	                }
	            }
	        }
	    }

	    private boolean isBorderCell(int row, int col, boolean includeDiagonals) {
	        Set<Integer> adjacentClusters = new HashSet<>();
	        int currentCluster = clusterGrid[row][col];
	        
	        // Check direct neighbors (up, down, left, right)
	        int[][] directNeighbors = {{-1,0}, {1,0}, {0,-1}, {0,1}};
	        for (int[] neighbor : directNeighbors) {
	            int newRow = row + neighbor[0];
	            int newCol = col + neighbor[1];
	            
	            if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols 
	                && clusterGrid[newRow][newCol] != -1
	                && clusterGrid[newRow][newCol] != currentCluster) {

	                adjacentClusters.add(clusterGrid[newRow][newCol]);
	            }
	        }
	        
	        // Check diagonal neighbors if requested
	        if (includeDiagonals) {
	            int[][] diagonalNeighbors = {{-1,-1}, {-1,1}, {1,-1}, {1,1}};
	            for (int[] neighbor : diagonalNeighbors) {
	                int newRow = row + neighbor[0];
	                int newCol = col + neighbor[1];
	                
	                if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols 
	                    && clusterGrid[newRow][newCol] != -1
	                    && clusterGrid[newRow][newCol] != currentCluster) {
	                    adjacentClusters.add(clusterGrid[newRow][newCol]);
	                }
	            }
	        }

	        if(!adjacentClusters.isEmpty()) {  // Changed condition: any different adjacent cluster makes it a border cell
	            Point p = new Point(row, col);
	            borderCells.add(p);
	            borderCellClusters.put(p, adjacentClusters);
	        }
	        return !adjacentClusters.isEmpty();  // Return true if there are any adjacent clusters
	    }

	    public void divideGridToDiagonalClusters(int numClusters) {
	        initializeClusters(0);
	        
	        // Calculate diagonal width to achieve desired number of clusters
	        int diagonalWidth = (int) Math.ceil((rows + cols) / (double) numClusters);
	        
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                // Assign cluster based on position along main diagonal
	                int clusterId = Math.min((i + j) / diagonalWidth, numClusters - 1);
	                clusterGrid[i][j] = clusterId;
	            }
	        }
	        
	        calculateBorderCells();
	    }

	    public void divideGridToCheckerboardClusters(int numClusters) {
	        initializeClusters(0);
	        
	        // Calculate checkerboard dimensions
	        int boardSize = (int) Math.ceil(Math.sqrt(numClusters));
	        int cellSize = Math.max(rows, cols) / boardSize;
	        
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                int blockRow = i / cellSize;
	                int blockCol = j / cellSize;
	                int clusterId = (blockRow * boardSize + blockCol) % numClusters;
	                clusterGrid[i][j] = clusterId;
	            }
	        }
	        
	        calculateBorderCells();
	    }

	    public void divideGridToFractalClusters(int numClusters) {
	        initializeClusters(0);
	        
	        // Use Hilbert curve for fractal-like partitioning
	        int order = (int) Math.ceil(Math.log(numClusters) / Math.log(4));
	        int totalPoints = 1 << (2 * order);
	        double pointsPerCluster = totalPoints / (double) numClusters;
	        
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                // Map 2D coordinates to Hilbert curve position
	                int hilbertPos = getHilbertCurvePosition(i * (1 << order) / rows, 
	                                                       j * (1 << order) / cols, 
	                                                       order);
	                int clusterId = Math.min((int)(hilbertPos / pointsPerCluster), numClusters - 1);
	                clusterGrid[i][j] = clusterId;
	            }
	        }
	        
	        calculateBorderCells();
	    }

	    public void divideGridToHoneycombClusters(int numClusters) {
	        initializeClusters(0);
	        
	        // Calculate hexagon size to achieve desired number of clusters
	        double hexSize = Math.sqrt((rows * cols) / (numClusters * 2.598)); // 2.598 is hex area ratio
	        
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                // Convert to hex coordinates
	                double q = (2.0/3 * j) / hexSize;
	                double r = (-1.0/3 * j + Math.sqrt(3)/3 * i) / hexSize;
	                
	                // Get nearest hex center
	                int qGrid = (int) Math.round(q);
	                int rGrid = (int) Math.round(r);
	                
	                // Convert hex coordinate to cluster ID
	                int clusterId = Math.abs((qGrid * 31 + rGrid * 37)) % numClusters;
	                clusterGrid[i][j] = clusterId;
	            }
	        }
	        
	        calculateBorderCells();
	    }

	    public void divideGridToWaveClusters(int numClusters) {
	        initializeClusters(0);
	        
	        // Use sine waves with different phases for interesting patterns
	        double frequency = 2 * Math.PI * 3 / Math.max(rows, cols);
	        
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                double wave1 = Math.sin(frequency * i);
	                double wave2 = Math.cos(frequency * j);
	                double combined = (wave1 + wave2 + 2) / 4; // Normalize to 0-1
	                
	                int clusterId = Math.min((int)(combined * numClusters), numClusters - 1);
	                clusterGrid[i][j] = clusterId;
	            }
	        }
	        
	        calculateBorderCells();
	    }

    public void divideGridToWeightedClusters(int numClusters, Integer[][] agentMask) {
        initializeClusters(numClusters);
        
        // Set the score grid from the agent mask
        setScoreGrid(agentMask);
        
        // Collect all agent cells with their scores
        List<AgentCell> agentCells = new ArrayList<>();
        int totalScore = 0;
        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (scoreGrid[i][j] > 0) {
                    AgentCell cell = new AgentCell(i, j, scoreGrid[i][j]);
                    agentCells.add(cell);
                    totalScore += scoreGrid[i][j];
                }
            }
        }
        
        if (agentCells.isEmpty()) {
            // If no agent cells, just use regular grid partitioning
            divideGridToGridClusters(numClusters);
            return;
        }
        
        // Calculate target score per cluster
        int baseScorePerCluster = totalScore / numClusters;
        int remainingScore = totalScore % numClusters;
        
        // Create seeds using spatial subdivision for better distribution
        List<AgentCell> seeds = createWeightedSeeds(agentCells, numClusters);
        
        // Initialize cluster assignments for seeds
        Map<Integer, Integer> clusterScores = new HashMap<>();
        for (int i = 0; i < numClusters; i++) {
            clusterScores.put(i, 0);
        }
        
        for (int i = 0; i < seeds.size(); i++) {
            AgentCell seed = seeds.get(i);
            clusterGrid[seed.row][seed.col] = i;
            clusterScores.put(i, clusterScores.get(i) + seed.score);
        }
        
        // Assign remaining agents using score-based load balancing
        List<AgentCell> unassignedCells = new ArrayList<>();
        for (AgentCell cell : agentCells) {
            if (clusterGrid[cell.row][cell.col] == -1) {
                unassignedCells.add(cell);
            }
        }
        
        // Sort by score (higher scores first) for better distribution
        unassignedCells.sort((a, b) -> Integer.compare(b.score, a.score));
        
        for (AgentCell cell : unassignedCells) {
            int bestCluster = findBestClusterForCell(cell, numClusters, clusterScores, 
                                                   baseScorePerCluster, remainingScore);
            
            if (bestCluster != -1) {
                clusterGrid[cell.row][cell.col] = bestCluster;
                clusterScores.put(bestCluster, clusterScores.get(bestCluster) + cell.score);
            }
        }
        
        // Refinement phase to improve spatial coherence for agent cells
        for (int iteration = 0; iteration < 50; iteration++) {
            boolean improved = false;
            Collections.shuffle(agentCells, random);
            
            for (AgentCell cell : agentCells) {
                if (clusterGrid[cell.row][cell.col] == -1) continue;
                
                int currentCluster = clusterGrid[cell.row][cell.col];
                Set<Integer> neighborClusters = getNeighborClusters(cell.row, cell.col);
                
                for (int newCluster : neighborClusters) {
                    if (newCluster == currentCluster) continue;
                    
                    double currentCoherence = calculateSpatialCoherence(cell.row, cell.col, currentCluster);
                    double newCoherence = calculateSpatialCoherence(cell.row, cell.col, newCluster);
                    
                    if (newCoherence > currentCoherence + 0.1) {
                        clusterGrid[cell.row][cell.col] = newCluster;
                        clusterScores.put(currentCluster, clusterScores.get(currentCluster) - cell.score);
                        clusterScores.put(newCluster, clusterScores.get(newCluster) + cell.score);
                        improved = true;
                        break;
                    }
                }
            }
            
            if (!improved) break;
        }
        
        // NEW: Assign all remaining cells (score 0) to nearest cluster
        assignEmptyCellsToNearestCluster();
        
        calculateBorderCells();
    }	    private List<AgentCell> createWeightedSeeds(List<AgentCell> agentCells, int numClusters) {
	        List<AgentCell> seeds = new ArrayList<>();
	        
	        if (numClusters <= agentCells.size()) {
	            // Use recursive spatial subdivision
	            List<int[]> regions = new ArrayList<>();
	            regions.add(new int[]{0, 0, rows, cols}); // [startRow, startCol, endRow, endCol]
	            
	            while (regions.size() < numClusters && !regions.isEmpty()) {
	                List<int[]> newRegions = new ArrayList<>();
	                for (int[] region : regions) {
	                    int startRow = region[0], startCol = region[1];
	                    int endRow = region[2], endCol = region[3];
	                    int midRow = (startRow + endRow) / 2;
	                    int midCol = (startCol + endCol) / 2;
	                    
	                    if (endRow - startRow > endCol - startCol) {
	                        newRegions.add(new int[]{startRow, startCol, midRow, endCol});
	                        newRegions.add(new int[]{midRow, startCol, endRow, endCol});
	                    } else {
	                        newRegions.add(new int[]{startRow, startCol, endRow, midCol});
	                        newRegions.add(new int[]{startRow, midCol, endRow, endCol});
	                    }
	                }
	                regions = newRegions;
	                if (regions.size() >= numClusters) break;
	            }
	            
	            // Find best agent in each region
	            for (int i = 0; i < numClusters && i < regions.size(); i++) {
	                int[] region = regions.get(i);
	                int centerRow = (region[0] + region[2]) / 2;
	                int centerCol = (region[1] + region[3]) / 2;
	                
	                AgentCell bestAgent = null;
	                double bestScore = -1;
	                
	                for (AgentCell agent : agentCells) {
	                    if (seeds.contains(agent)) continue;
	                    
	                    boolean inRegion = agent.row >= region[0] && agent.row < region[2] && 
	                                     agent.col >= region[1] && agent.col < region[3];
	                    double dist = Math.sqrt(Math.pow(agent.row - centerRow, 2) + 
	                                          Math.pow(agent.col - centerCol, 2));
	                    
	                    // Score combines agent score with spatial preference
	                    double score = agent.score * 2.0 - dist * 0.1;
	                    if (inRegion) score += 10.0; // Strong preference for agents in region
	                    
	                    if (score > bestScore) {
	                        bestScore = score;
	                        bestAgent = agent;
	                    }
	                }
	                
	                if (bestAgent != null) {
	                    seeds.add(bestAgent);
	                }
	            }
	        }
	        
	        return seeds;
	    }
	    
	    private int findBestClusterForCell(AgentCell cell, int numClusters, 
	                                     Map<Integer, Integer> clusterScores,
	                                     int baseScorePerCluster, int remainingScore) {
	        int bestCluster = -1;
	        double bestScore = Double.MAX_VALUE;
	        
	        for (int cluster = 0; cluster < numClusters; cluster++) {
	            int currentScore = clusterScores.get(cluster);
	            int targetScore = baseScorePerCluster + (cluster < remainingScore ? 1 : 0);
	            
	            if (currentScore >= targetScore + cell.score) continue; // Would overflow
	            
	            // Calculate spatial distance to cluster
	            double minDistToCluster = Double.MAX_VALUE;
	            for (int i = 0; i < rows; i++) {
	                for (int j = 0; j < cols; j++) {
	                    if (clusterGrid[i][j] == cluster) {
	                        double dist = Math.sqrt(Math.pow(cell.row - i, 2) + Math.pow(cell.col - j, 2));
	                        minDistToCluster = Math.min(minDistToCluster, dist);
	                    }
	                }
	            }
	            
	            // Score combines distance and load balancing
	            double loadFactor = (double)currentScore / targetScore;
	            double score = minDistToCluster * (1.0 + loadFactor);
	            
	            if (score < bestScore) {
	                bestScore = score;
	                bestCluster = cluster;
	            }
	        }
	        
	        return bestCluster;
	    }
	    
	    private Set<Integer> getNeighborClusters(int i, int j) {
	        Set<Integer> neighbors = new HashSet<>();
	        int[][] dirs = {{-1,0}, {1,0}, {0,-1}, {0,1}};
	        for (int[] dir : dirs) {
	            int ni = i + dir[0];
	            int nj = j + dir[1];
	            if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && 
	                clusterGrid[ni][nj] >= 0) {
	                neighbors.add(clusterGrid[ni][nj]);
	            }
	        }
	        return neighbors;
	    }
	    
    private double calculateSpatialCoherence(int i, int j, int targetCluster) {
        int sameClusterNeighbors = 0;
        int totalNeighbors = 0;
        
        for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
                if (di == 0 && dj == 0) continue;
                
                int ni = i + di;
                int nj = j + dj;
                if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && 
                    clusterGrid[ni][nj] >= 0) {
                    totalNeighbors++;
                    if (clusterGrid[ni][nj] == targetCluster) {
                        sameClusterNeighbors++;
                    }
                }
            }
        }
        
        return totalNeighbors > 0 ? (double)sameClusterNeighbors / totalNeighbors : 0;
    }
    
    /**
     * Assigns all empty cells (score = 0) to the nearest cluster containing agent cells
     */
    private void assignEmptyCellsToNearestCluster() {
        // Use a queue for breadth-first search to assign empty cells
        Queue<Point> queue = new LinkedList<>();
        
        // Add all cells that are already assigned to clusters to the queue
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (clusterGrid[i][j] != -1) {
                    queue.offer(new Point(i, j));
                }
            }
        }
        
        // BFS to assign unassigned cells to nearest cluster
        int[][] directions = {{-1,0}, {1,0}, {0,-1}, {0,1}};
        
        while (!queue.isEmpty()) {
            Point current = queue.poll();
            int currentCluster = clusterGrid[current.row][current.col];
            
            // Check all neighbors
            for (int[] dir : directions) {
                int newRow = current.row + dir[0];
                int newCol = current.col + dir[1];
                
                if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols 
                    && clusterGrid[newRow][newCol] == -1) {
                    
                    // Assign this cell to the same cluster as current
                    clusterGrid[newRow][newCol] = currentCluster;
                    queue.offer(new Point(newRow, newCol));
                }
            }
        }
        
        // Fallback: if any cells are still unassigned, assign them to cluster 0
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (clusterGrid[i][j] == -1) {
                    clusterGrid[i][j] = 0;
                }
            }
        }
    }	    private int getHilbertCurvePosition(int x, int y, int order) {
	        int position = 0;
	        for (int i = 0; i < order; i++) {
	            int xi = (x >> i) & 1;
	            int yi = (y >> i) & 1;
	            position += ((3 * xi) ^ yi) << (2 * i);
	        }
	        return position;
	    }

	    private void calculateBorderCells() {
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }


	    public Set<Integer> getBorderCellClusters(Point p) {
	        return borderCellClusters.getOrDefault(p, new HashSet<>());
	    }

	    private int getCellId(int row, int col) {
	        return row * cols + col;
	    }

    public List<List<Integer>> getClusteringAsLists() {
        if (clusterGrid == null) {
            return new ArrayList<>();
        }

        Map<Integer, List<Integer>> clusterMap = new HashMap<>();
        
        // Include ALL cells in clusters, regardless of their score
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (clusterGrid[i][j] != -1) {
                    int clusterId = clusterGrid[i][j];
                    int cellId = getCellId(i, j);
                    clusterMap.computeIfAbsent(clusterId, k -> new ArrayList<>()).add(cellId);
                }
            }
        }
        
        System.out.println("clusterMap " + clusterMap);
        System.out.println("ArrayList<>(clusterMap.values() " + new ArrayList<>(clusterMap.values()));
        // Convert map to list of lists
        return new ArrayList<>(clusterMap.values());
    }	    public void divideGridToGridClusters(int numClusters) {
	        initializeClusters(0);

	        // Find the best grid dimensions that are close to square
	        int numRows = (int) Math.sqrt(numClusters);
	        while (numClusters % numRows != 0 && numRows > 1) {
	            numRows--;
	        }
	        int numCols = numClusters / numRows;

	        // Calculate the size of each grid cell
	        int baseHeight = rows / numRows;
	        int baseWidth = cols / numCols;
	        
	        // Calculate remainders
	        int remainderRows = rows % numRows;
	        int remainderCols = cols % numCols;

	        // Create arrays to store the actual size of each row and column
	        int[] rowHeights = new int[numRows];
	        int[] colWidths = new int[numCols];

	        // Distribute the remainder pixels evenly
	        for (int i = 0; i < numRows; i++) {
	            rowHeights[i] = baseHeight + (i < remainderRows ? 1 : 0);
	        }
	        for (int i = 0; i < numCols; i++) {
	            colWidths[i] = baseWidth + (i < remainderCols ? 1 : 0);
	        }

	        // Assign clusters based on position
	        int currentRow = 0;
	        int currentY = 0;
	        
	        for (int i = 0; i < rows; i++) {
	            if (i >= currentY + rowHeights[currentRow]) {
	                currentY += rowHeights[currentRow];
	                currentRow++;
	            }
	            
	            int currentCol = 0;
	            int currentX = 0;
	            
	            for (int j = 0; j < cols; j++) {
	                if (j >= currentX + colWidths[currentCol]) {
	                    currentX += colWidths[currentCol];
	                    currentCol++;
	                }
	                
	                clusterGrid[i][j] = currentRow * numCols + currentCol;
	            }
	        }

	        // Calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }

	    public void divideGridToHorizontalStrips(int numClusters) {
	        initializeClusters(0);

	        // Calculate the height of each strip
	        int stripHeight = (int) Math.ceil(rows / (double) numClusters);

	        // Assign each cell to its strip cluster
	        for (int i = 0; i < rows; i++) {
	            int clusterId = Math.min(i / stripHeight, numClusters - 1);
	            for (int j = 0; j < cols; j++) {
	                clusterGrid[i][j] = clusterId;
	            }
	        }

	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }

	    public void divideGridToVerticalStrips(int numClusters) {
	        initializeClusters(0);

	        // Calculate the width of each strip
	        int stripWidth = (int) Math.ceil(cols / (double) numClusters);

	        // Assign each cell to its strip cluster
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                int clusterId = Math.min(j / stripWidth, numClusters - 1);
	                clusterGrid[i][j] = clusterId;
	            }
	        }

	        // Calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }

	    private void recursiveSplit(int startRow, int startCol, int endRow, int endCol, 
	                              int currentCluster, int remainingClusters, boolean splitVertical) {
	        // Base case: if we only need one cluster, assign all cells to it
	        if (remainingClusters == 1) {
	            for (int i = startRow; i < endRow; i++) {
	                for (int j = startCol; j < endCol; j++) {
	                    clusterGrid[i][j] = currentCluster;
	                }
	            }
	            return;
	        }

	        // Calculate how to split remaining clusters
	        int firstPartClusters = remainingClusters / 2;
	        int secondPartClusters = remainingClusters - firstPartClusters;

	        if (splitVertical) {
	            // Split vertically
	            int midCol = startCol + (endCol - startCol) * firstPartClusters / remainingClusters;
	            recursiveSplit(startRow, startCol, endRow, midCol, currentCluster, firstPartClusters, false);
	            recursiveSplit(startRow, midCol, endRow, endCol, currentCluster + firstPartClusters, secondPartClusters, false);
	        } else {
	            // Split horizontally
	            int midRow = startRow + (endRow - startRow) * firstPartClusters / remainingClusters;
	            recursiveSplit(startRow, startCol, midRow, endCol, currentCluster, firstPartClusters, true);
	            recursiveSplit(midRow, startCol, endRow, endCol, currentCluster + firstPartClusters, secondPartClusters, true);
	        }
	    }

	    public void divideGridToBSPClusters(int numClusters) {
	        initializeClusters(0);

	        // Start recursive splitting
	        recursiveSplit(0, 0, rows, cols, 0, numClusters, true);

	        // Calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }

	    public void divideGridToCircularClusters(int numClusters) {
	        initializeClusters(0);
	        
	        // Calculate center point
	        int centerX = cols / 2;
	        int centerY = rows / 2;
	        
	        // Calculate max radius and radius step
	        double maxRadius = Math.sqrt(Math.pow(Math.max(rows, cols), 2) / 2);
	        double radiusStep = maxRadius / numClusters;
	        
	        // Assign clusters based on distance from center
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                double distance = Math.sqrt(Math.pow(i - centerY, 2) + Math.pow(j - centerX, 2));
	                int clusterId = Math.min((int)(distance / radiusStep), numClusters - 1);
	                clusterGrid[i][j] = clusterId;
	            }
	        }

	        // Calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }

	    public void divideGridToSpiralClusters(int numClusters) {
	        initializeClusters(0);
	        
	        int currentCluster = 0;
	        int left = 0, right = cols-1, top = 0, bottom = rows-1;
	        int cellsPerCluster = (rows * cols) / numClusters;
	        int cellsInCurrentCluster = 0;
	        
	        while (left <= right && top <= bottom) {
	            // Move right
	            for (int j = left; j <= right && currentCluster < numClusters; j++) {
	                clusterGrid[top][j] = currentCluster;
	                cellsInCurrentCluster++;
	                if (cellsInCurrentCluster >= cellsPerCluster && currentCluster < numClusters - 1) {
	                    currentCluster++;
	                    cellsInCurrentCluster = 0;
	                }
	            }
	            top++;
	            
	            // Move down
	            for (int i = top; i <= bottom && currentCluster < numClusters; i++) {
	                clusterGrid[i][right] = currentCluster;
	                cellsInCurrentCluster++;
	                if (cellsInCurrentCluster >= cellsPerCluster && currentCluster < numClusters - 1) {
	                    currentCluster++;
	                    cellsInCurrentCluster = 0;
	                }
	            }
	            right--;
	            
	            // Move left
	            for (int j = right; j >= left && currentCluster < numClusters; j--) {
	                clusterGrid[bottom][j] = currentCluster;
	                cellsInCurrentCluster++;
	                if (cellsInCurrentCluster >= cellsPerCluster && currentCluster < numClusters - 1) {
	                    currentCluster++;
	                    cellsInCurrentCluster = 0;
	                }
	            }
	            bottom--;
	            
	            // Move up
	            for (int i = bottom; i >= top && currentCluster < numClusters; i--) {
	                clusterGrid[i][left] = currentCluster;
	                cellsInCurrentCluster++;
	                if (cellsInCurrentCluster >= cellsPerCluster && currentCluster < numClusters - 1) {
	                    currentCluster++;
	                    cellsInCurrentCluster = 0;
	                }
	            }
	            left++;
	        }

	        // Fill any remaining cells with the last cluster
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (clusterGrid[i][j] == -1) {
	                    clusterGrid[i][j] = numClusters - 1;
	                }
	            }
	        }

	        // Calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }

	    private int findNearestNonEmptyCluster(int row, int col) {
	        int radius = 1;
	        while (radius < Math.max(rows, cols)) {
	            for (int i = -radius; i <= radius; i++) {
	                for (int j = -radius; j <= radius; j++) {
	                    if (Math.abs(i) == radius || Math.abs(j) == radius) {
	                        int ni = row + i;
	                        int nj = col + j;
	                        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && clusterGrid[ni][nj] != -1) {
	                            return clusterGrid[ni][nj];
	                        }
	                    }
	                }
	            }
	            radius++;
	        }
	        return 0;
	    }

	    public void divideGridToVoronoiClusters(int numClusters) {
	        initializeClusters(0);
	        
	        int targetSize = (rows * cols) / numClusters;
	        Map<Integer, Integer> clusterSizes = new HashMap<>();
	        PriorityQueue<Point> growthFrontier = new PriorityQueue<>((a, b) -> {
	            int clusterA = clusterGrid[a.row][a.col];
	            int clusterB = clusterGrid[b.row][b.col];
	            return Integer.compare(clusterSizes.get(clusterA), clusterSizes.get(clusterB));
	        });

	        // Place seeds in a roughly uniform distribution
	        double seedSpacing = Math.sqrt((rows * cols) / (double) numClusters);
	        for (int k = 0; k < numClusters; k++) {
	            int row = (int) ((k / Math.sqrt(numClusters)) * seedSpacing + random.nextDouble() * seedSpacing);
	            int col = (int) ((k % Math.sqrt(numClusters)) * seedSpacing + random.nextDouble() * seedSpacing);
	            
	            row = Math.min(Math.max(row, 0), rows - 1);
	            col = Math.min(Math.max(col, 0), cols - 1);
	            
	            clusterGrid[row][col] = k;
	            clusterSizes.put(k, 1);
	            growthFrontier.add(new Point(row, col));
	        }

	        // Grow regions while maintaining size balance
	        int[][] directions = {{-1,0}, {1,0}, {0,-1}, {0,1}};
	        boolean[][] visited = new boolean[rows][cols];
	        
	        while (!growthFrontier.isEmpty()) {
	            Point current = growthFrontier.poll();
	            int currentCluster = clusterGrid[current.row][current.col];
	            
	            // Skip if this cluster is already at target size
	            if (clusterSizes.get(currentCluster) >= targetSize) {
	                continue;
	            }

	            // Try to grow in each direction
	            for (int[] dir : directions) {
	                int newRow = current.row + dir[0];
	                int newCol = current.col + dir[1];
	                
	                if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols 
	                    && !visited[newRow][newCol] 
	                    && clusterGrid[newRow][newCol] == -1) {
	                    
	                    clusterGrid[newRow][newCol] = currentCluster;
	                    clusterSizes.merge(currentCluster, 1, Integer::sum);
	                    visited[newRow][newCol] = true;
	                    growthFrontier.add(new Point(newRow, newCol));
	                }
	            }
	        }

	        // Assign any remaining unassigned cells to nearest cluster
	        boolean changed;
	        do {
	            changed = false;
	            for (int i = 0; i < rows; i++) {
	                for (int j = 0; j < cols; j++) {
	                    if (clusterGrid[i][j] == -1) {
	                        int nearestCluster = findNearestNonEmptyCluster(i, j);
	                        clusterGrid[i][j] = nearestCluster;
	                        changed = true;
	                    }
	                }
	            }
	        } while (changed);

	        // Calculate border cells
	        for (int i = 0; i < rows; i++) {
	            for (int j = 0; j < cols; j++) {
	                if (isBorderCell(i, j, false)) {
	                    Point p = new Point(i, j);
	                    adjCell.add(p);
	                }
	            }
	        }
	    }
	}
	
	/**
     * Grid partitioning using Voronoi technique
     * 
     * @param rows : number of rows in the grid
     * @param columns : number of columns in the grid
     * @param cluster_number : number of clusters required
     * @param number_neighbors : number of neighbors to consider
     * 
     * @return List<List<Integer>> clustered graph
     */
    @SuppressWarnings("rawtypes")
    @operator(
            value = "grid_voronoi_partitioning",
            type = IType.LIST,
            category = { IOperatorCategory.GRID },
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Voronoi technique")
    public static List<List<Integer>> grid_voronoi_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToVoronoiClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Grid technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_grid_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Grid technique")
    public static List<List<Integer>> grid_grid_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToGridClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Horizontal Strips technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_horizontal_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Horizontal Strips technique")
    public static List<List<Integer>> grid_horizontal_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToHorizontalStrips(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Vertical Strips technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_vertical_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Vertical Strips technique")
    public static List<List<Integer>> grid_vertical_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToVerticalStrips(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using BSP technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_bsp_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using BSP technique")
    public static List<List<Integer>> grid_bsp_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToBSPClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Circular technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_circular_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Circular technique")
    public static List<List<Integer>> grid_circular_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToCircularClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Spiral technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_spiral_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Spiral technique")
    public static List<List<Integer>> grid_spiral_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToSpiralClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Diagonal technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_diagonal_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Diagonal technique")
    public static List<List<Integer>> grid_diagonal_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToDiagonalClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Checkerboard technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_checkerboard_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Checkerboard technique")
    public static List<List<Integer>> grid_checkerboard_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToCheckerboardClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Fractal technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_fractal_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Fractal technique")
    public static List<List<Integer>> grid_fractal_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToFractalClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Honeycomb technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_honeycomb_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Honeycomb technique")
    public static List<List<Integer>> grid_honeycomb_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToHoneycombClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Grid partitioning using Wave technique
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_wave_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning using Wave technique")
    public static List<List<Integer>> grid_wave_partitioning(IScope scope, int rows, int columns, int cluster_number, int number_neighbors) {
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        try {
            grid.divideGridToWaveClusters(cluster_number);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * Helper method to convert ArrayList<Integer> (1D list) to Integer[][] (2D grid)
     * Assumes the 1D list represents a flattened 2D grid in row-major order
     */
    private static Integer[][] convertToGrid(ArrayList<Integer> agentMask, int rows, int cols) {
        if (agentMask == null || agentMask.isEmpty()) {
            return new Integer[rows][cols]; // Return empty grid with specified dimensions
        }
        
        Integer[][] grid = new Integer[rows][cols];
        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                int index = i * cols + j; // Convert 2D coordinates to 1D index
                grid[i][j] = (index < agentMask.size()) ? agentMask.get(index) : 0;
            }
        }
        
        return grid;
    }

    /**
     * Grid partitioning technique with scoring system
     */
    @SuppressWarnings("rawtypes")
    @operator(value = "grid_score_partitioning", 
    		type = IType.LIST, 
            category = { IOperatorCategory.GRID }, 
            concept = { IConcept.GRID })
    @doc(value = "grid partitioning technique with scoring system. The agentMask parameter should be a 1D list representing a flattened 2D grid where 0=no agent, 1=bacteria, 2=nematodes")
    public static List<List<Integer>> grid_score_partitioning(IScope scope,int rows, int columns, ArrayList<Integer> agentMask, int cluster_number) {
        if (agentMask == null || agentMask.isEmpty()) {
            return new ArrayList<>();
        }
        
        GRID grid = new GRID(rows, columns);
        random = scope.getRandom().getGenerator();
        
        try {
            // Convert 1D ArrayList to 2D array format for internal processing
            Integer[][] maskArray = convertToGrid(agentMask, rows, columns);
            grid.divideGridToWeightedClusters(cluster_number, maskArray);
            return grid.getClusteringAsLists();
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    /**
     * 
     * Grid partitionning of a grid in cluster_number cluster(s)
     * 
     * 
     * @param edges : edges of the graph to clsuter
     * @param cluster_number : number of cluster requireds
     * 
     * @return List<List<Integer>> clustered graph
     */
    @SuppressWarnings("rawtypes")
	@operator (
			value = "grid_KMEAN_partitionning",
			type = IType.LIST,
			category = { IOperatorCategory.GRID },
			concept = { IConcept.GRID})
	@doc (value = "grid partitionning of a grid")
    public static List<List<Integer>> grid_KMEAN_partitionning(IScope scope,int rows, int columns, int cluster_number, int number_neighbors) 
    {
    	GRID grid = new GRID(rows, columns);
    	System.out.println("scope " + scope);
    	
    	random = scope.getRandom().getGenerator();
        // Choose partitioning method based on analysis results
        try {
            grid.divideGridToKMEANClusters(cluster_number);
            List<List<Integer>> clusters = grid.getClusteringAsLists();
            return clusters;
        } catch (Exception e) {
            e.printStackTrace();
            return new ArrayList<>();
        }
    }
}
