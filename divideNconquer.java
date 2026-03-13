import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;
import java.util.Timer;

public class TracksGame extends JFrame {

    private int gridSize = 8;
    private static final int CELL_SIZE = 50;

    private GamePanel gamePanel;
    private SidePanel sidePanel;
    private GameMode currentMode = GameMode.PATTERN1;

    enum GameMode { PATTERN1, PATTERN2, PATTERN3, PATTERN4, PATTERN5, RANDOM }

    public TracksGame() {
        setTitle("Tracks Lab - Quadtree Divide & Conquer");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        gamePanel = new GamePanel();
        sidePanel = new SidePanel(gamePanel);

        add(gamePanel, BorderLayout.CENTER);
        add(sidePanel, BorderLayout.EAST);

        pack();
        setLocationRelativeTo(null);
        setResizable(false);
    }

    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception ignored) {}
        SwingUtilities.invokeLater(() -> new TracksGame().setVisible(true));
    }

    /**
     * Cell class representing each grid position
     */
    static class Cell {
        int r, c;
        boolean fromPlayer = false;
        boolean fromSolver = false;
        boolean isCrossed = false;
        boolean isHint = false;
        
        Cell playerParent = null;
        Cell solverParent = null;

        public Cell(int r, int c) {
            this.r = r;
            this.c = c;
        }

        boolean isUsed() {
            return fromPlayer || fromSolver;
        }

        boolean isOwnedByPlayer() {
            return fromPlayer;
        }

        boolean isOwnedBySolver() {
            return fromSolver;
        }

        @Override
        public String toString() {
            return "(" + r + "," + c + ")";
        }
    }

    /**
     * Quadrant result class for D&C analysis
     */
    class QuadrantResult {
        int r1, r2, c1, c2;
        int score;
        int[] rowUsed;
        int[] colUsed;
        boolean[] topBoundary;
        boolean[] bottomBoundary;
        boolean[] leftBoundary;
        boolean[] rightBoundary;
        int usedCells;
        int exactRowMatches;
        int exactColMatches;
        int boundaryViolations;
        int continuityScore;

        QuadrantResult(int r1, int r2, int c1, int c2) {
            this.r1 = r1;
            this.r2 = r2;
            this.c1 = c1;
            this.c2 = c2;
            int rows = r2 - r1 + 1;
            int cols = c2 - c1 + 1;
            rowUsed = new int[rows];
            colUsed = new int[cols];
            topBoundary = new boolean[cols];
            bottomBoundary = new boolean[cols];
            leftBoundary = new boolean[rows];
            rightBoundary = new boolean[rows];
        }
    }

    /**
     * Main game panel with all game logic
     */
    class GamePanel extends JPanel {
        private Cell[][] grid;
        private int[] rowClues, colClues;
        private Cell startNode, endNode;
        private Cell playerHead, solverHead;
        private boolean isPlayerTurn = true;
        private boolean gameOver = false;
        private boolean computerThinking = false;
        private boolean autoSolving = false;

        public GamePanel() {
            setBackground(Color.WHITE);
            setPreferredSize(new Dimension(gridSize * CELL_SIZE + 100, gridSize * CELL_SIZE + 100));
            initGame(GameMode.PATTERN1);

            addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent e) {
                    if (gameOver || !isPlayerTurn || computerThinking || autoSolving) return;

                    int startX = 50, startY = 50;
                    int c = (e.getX() - startX) / CELL_SIZE;
                    int r = (e.getY() - startY) / CELL_SIZE;

                    if (r < 0 || c < 0 || r >= gridSize || c >= gridSize) return;

                    if (SwingUtilities.isRightMouseButton(e)) {
                        // Right click to place X
                        Cell clicked = grid[r][c];
                        if (!clicked.isUsed() && clicked != startNode && clicked != endNode) {
                            clicked.isCrossed = !clicked.isCrossed;
                            repaint();
                        }
                        return;
                    }

                    if (SwingUtilities.isLeftMouseButton(e) && isValidPlayerMove(r, c)) {
                        // Player makes a move
                        clearHints();
                        makePlayerMove(r, c);
                        repaint();
                        checkGameStatus();

                        if (!gameOver) {
                            // Computer's turn
                            computerThinking = true;
                            isPlayerTurn = false;
                            sidePanel.updateStatus("Computer thinking (Quadrant Analysis)...");
                            
                            // Add small delay for better UX
                            new Timer().schedule(new TimerTask() {
                                @Override
                                public void run() {
                                    computerTurn();
                                }
                            }, 150);
                        }
                    }
                }
            });
        }

        /**
         * Initialize game with selected pattern
         */
        public void initGame(GameMode mode) {
            currentMode = mode;

            switch (mode) {
                case PATTERN1: gridSize = 8; break;
                case PATTERN2: gridSize = 9; break;
                case PATTERN3: gridSize = 10; break;
                case PATTERN4: gridSize = 8; break;
                case PATTERN5: gridSize = 10; break;
                case RANDOM: gridSize = 8; break;
            }

            setPreferredSize(new Dimension(gridSize * CELL_SIZE + 100, gridSize * CELL_SIZE + 100));
            grid = new Cell[gridSize][gridSize];
            
            for (int r = 0; r < gridSize; r++) {
                for (int c = 0; c < gridSize; c++) {
                    grid[r][c] = new Cell(r, c);
                }
            }

            if (mode == GameMode.PATTERN1) {
                generatePattern1();
            } else {
                long seed;
                switch (mode) {
                    case PATTERN2: seed = 12345L; break;
                    case PATTERN3: seed = 67890L; break;
                    case PATTERN4: seed = 11111L; break;
                    case PATTERN5: seed = 22222L; break;
                    default: seed = System.currentTimeMillis(); break;
                }
                generateComplexPuzzle(seed);
            }

            playerHead = startNode;
            solverHead = endNode;
            isPlayerTurn = true;
            gameOver = false;
            computerThinking = false;
            autoSolving = false;

            TracksGame.this.pack();
            repaint();

            if (sidePanel != null) {
                String modeName = (mode == GameMode.RANDOM) ? "Random" : "Pattern " + mode.name().charAt(7);
                sidePanel.updateStatus("Player's Turn (A → ...)");
                sidePanel.updateAnalysis(
                    "<html><b>Puzzle:</b> " + modeName +
                    "<br/>Size: " + gridSize + "x" + gridSize +
                    "<br/>Required Cells: " + sumArray(rowClues) +
                    "<br/><b>Algorithm:</b> Quadtree D&C<br/><b>Computer:</b> Continues from B</html>"
                );
            }
        }

        /**
         * Generate Pattern 1 (fixed path)
         */
        private void generatePattern1() {
            List<Cell> path = new ArrayList<>();
            int[][] coords = {
                {1,1}, {1,2}, {1,3}, {1,4}, {2,4}, {3,4},
                {3,3}, {3,2}, {3,1}, {4,1}, {5,1},
                {5,2}, {5,3}, {5,4}, {5,5}, {5,6},
                {4,6}, {3,6}, {2,6}, {1,6}
            };

            for (int[] p : coords) path.add(grid[p[0]][p[1]]);

            startNode = path.get(0);
            endNode = path.get(path.size() - 1);
            startNode.fromPlayer = true;
            endNode.fromSolver = true;

            calculateCluesFromPath(path);
        }

        /**
         * Generate complex puzzle with given seed
         */
        private void generateComplexPuzzle(long seed) {
            Random rand = new Random(seed);
            int targetLength = (int) ((gridSize * gridSize) * 0.55);

            while (true) {
                for (int r = 0; r < gridSize; r++) {
                    for (int c = 0; c < gridSize; c++) {
                        grid[r][c] = new Cell(r, c);
                    }
                }

                int startR = rand.nextInt(gridSize - 2) + 1;
                int startC = rand.nextInt(gridSize - 2) + 1;
                Cell start = grid[startR][startC];

                boolean[][] visited = new boolean[gridSize][gridSize];
                visited[startR][startC] = true;

                List<Cell> bestPath = new ArrayList<>();
                List<Cell> currentPath = new ArrayList<>();
                currentPath.add(start);

                generatePathDFS(start, visited, currentPath, bestPath, targetLength, rand);

                if (bestPath.size() >= targetLength * 0.8) {
                    startNode = bestPath.get(0);
                    endNode = bestPath.get(bestPath.size() - 1);
                    startNode.fromPlayer = true;
                    endNode.fromSolver = true;
                    calculateCluesFromPath(bestPath);
                    break;
                }
            }
        }

        /**
         * DFS for path generation
         */
        private void generatePathDFS(Cell current, boolean[][] visited, List<Cell> currentPath,
                                     List<Cell> bestPath, int target, Random rand) {
            if (currentPath.size() > bestPath.size()) {
                bestPath.clear();
                bestPath.addAll(currentPath);
            }
            if (bestPath.size() >= target) return;

            List<Cell> neighbors = getNeighbors(current);
            Collections.shuffle(neighbors, rand);

            for (Cell n : neighbors) {
                if (!visited[n.r][n.c] && countPathNeighbors(n, visited) <= 1) {
                    visited[n.r][n.c] = true;
                    currentPath.add(n);
                    generatePathDFS(n, visited, currentPath, bestPath, target, rand);
                    if (bestPath.size() >= target) return;
                    currentPath.remove(currentPath.size() - 1);
                    visited[n.r][n.c] = false;
                }
            }
        }

        /**
         * Count neighbors in path
         */
        private int countPathNeighbors(Cell c, boolean[][] visited) {
            int count = 0;
            for (Cell n : getNeighbors(c)) {
                if (visited[n.r][n.c]) count++;
            }
            return count;
        }

        /**
         * Calculate row/column clues from path
         */
        private void calculateCluesFromPath(List<Cell> path) {
            rowClues = new int[gridSize];
            colClues = new int[gridSize];
            for (Cell c : path) {
                rowClues[c.r]++;
                colClues[c.c]++;
            }
        }

        /**
         * Reset current puzzle
         */
        public void resetCurrentPuzzle() {
            initGame(currentMode);
        }

        /**
         * New random puzzle
         */
        public void newRandomPuzzle() {
            initGame(GameMode.RANDOM);
        }

        /**
         * Load pattern by number
         */
        public void loadPattern(int patternNum) {
            switch (patternNum) {
                case 1: initGame(GameMode.PATTERN1); break;
                case 2: initGame(GameMode.PATTERN2); break;
                case 3: initGame(GameMode.PATTERN3); break;
                case 4: initGame(GameMode.PATTERN4); break;
                case 5: initGame(GameMode.PATTERN5); break;
            }
        }

        /**
         * Get valid neighboring cells
         */
        public List<Cell> getNeighbors(Cell c) {
            List<Cell> neighbors = new ArrayList<>();
            int[] dr = {-1, 1, 0, 0};
            int[] dc = {0, 0, -1, 1};

            for (int i = 0; i < 4; i++) {
                int nr = c.r + dr[i];
                int nc = c.c + dc[i];
                if (nr >= 0 && nr < gridSize && nc >= 0 && nc < gridSize) {
                    neighbors.add(grid[nr][nc]);
                }
            }
            return neighbors;
        }

        /**
         * Manhattan distance
         */
        private int manhattan(Cell a, Cell b) {
            return Math.abs(a.r - b.r) + Math.abs(a.c - b.c);
        }

        /**
         * Sum array elements
         */
        private int sumArray(int[] arr) {
            int sum = 0;
            for (int v : arr) sum += v;
            return sum;
        }

        /**
         * Count used cells in a row
         */
        private int countUsedInRow(int r) {
            int count = 0;
            for (int c = 0; c < gridSize; c++) {
                if (grid[r][c].isUsed()) count++;
            }
            return count;
        }

        /**
         * Count used cells in a column
         */
        private int countUsedInCol(int c) {
            int count = 0;
            for (int r = 0; r < gridSize; r++) {
                if (grid[r][c].isUsed()) count++;
            }
            return count;
        }

        /**
         * Check if a player move is valid
         */
        public boolean isValidPlayerMove(int r, int c) {
            if (r < 0 || c < 0 || r >= gridSize || c >= gridSize) return false;

            Cell target = grid[r][c];
            if (target.isUsed() || target.isCrossed) return false;
            if (Math.abs(r - playerHead.r) + Math.abs(c - playerHead.c) != 1) return false;
            if (countUsedInRow(r) >= rowClues[r]) return false;
            if (countUsedInCol(c) >= colClues[c]) return false;

            return true;
        }

        /**
         * Check if a solver move is valid
         */
        public boolean isValidSolverMove(Cell target) {
            if (target == null) return false;
            if (target.isUsed() || target.isCrossed) return false;
            if (Math.abs(target.r - solverHead.r) + Math.abs(target.c - solverHead.c) != 1) return false;
            if (countUsedInRow(target.r) >= rowClues[target.r]) return false;
            if (countUsedInCol(target.c) >= colClues[target.c]) return false;

            return true;
        }

        /**
         * Make a player move
         */
        private void makePlayerMove(int r, int c) {
            Cell next = grid[r][c];
            next.fromPlayer = true;
            next.playerParent = playerHead;
            playerHead = next;
        }

        /**
         * Make a solver move
         */
        private void makeSolverMove(Cell next) {
            next.fromSolver = true;
            next.solverParent = solverHead;
            solverHead = next;
        }

        /**
         * Undo player move
         */
        private void undoPlayerMove(Cell cell, Cell oldHead) {
            cell.fromPlayer = false;
            cell.playerParent = null;
            playerHead = oldHead;
        }

        /**
         * Undo solver move
         */
        private void undoSolverMove(Cell cell, Cell oldHead) {
            cell.fromSolver = false;
            cell.solverParent = null;
            solverHead = oldHead;
        }

        /**
         * Check if heads are adjacent
         */
        private boolean headsAdjacent() {
            return Math.abs(playerHead.r - solverHead.r) + Math.abs(playerHead.c - solverHead.c) == 1;
        }

        /**
         * Count total used cells
         */
        private int totalCellsUsed() {
            int used = 0;
            for (int r = 0; r < gridSize; r++) {
                for (int c = 0; c < gridSize; c++) {
                    if (grid[r][c].isUsed()) used++;
                }
            }
            return used;
        }

        /**
         * QUADTREE DIVIDE & CONQUER - Core Algorithm
         * 
         * Recursively divides the board into quadrants and analyzes boundary conditions
         */
        private QuadrantResult analyzeQuadrant(int r1, int r2, int c1, int c2) {
            if (r1 > r2 || c1 > c2) return null;

            int height = r2 - r1 + 1;
            int width = c2 - c1 + 1;

            // Base case: small quadrant (2x2 or smaller)
            if (height <= 2 || width <= 2) {
                return analyzeBaseQuadrant(r1, r2, c1, c2);
            }

            // Divide into 4 quadrants
            int rm = (r1 + r2) / 2;
            int cm = (c1 + c2) / 2;

            QuadrantResult q1 = analyzeQuadrant(r1, rm, c1, cm);       // Top-left
            QuadrantResult q2 = analyzeQuadrant(r1, rm, cm + 1, c2);   // Top-right
            QuadrantResult q3 = analyzeQuadrant(rm + 1, r2, c1, cm);   // Bottom-left
            QuadrantResult q4 = analyzeQuadrant(rm + 1, r2, cm + 1, c2); // Bottom-right

            // Combine step - synchronize boundaries
            return combineQuadrants(r1, r2, c1, c2, q1, q2, q3, q4);
        }

        /**
         * Analyze base quadrant (2x2 or smaller)
         */
        private QuadrantResult analyzeBaseQuadrant(int r1, int r2, int c1, int c2) {
            QuadrantResult result = new QuadrantResult(r1, r2, c1, c2);
            int score = 0;

            // Count used cells and set boundary flags
            for (int r = r1; r <= r2; r++) {
                int localR = r - r1;
                for (int c = c1; c <= c2; c++) {
                    int localC = c - c1;
                    Cell cell = grid[r][c];
                    
                    if (cell.isUsed()) {
                        result.rowUsed[localR]++;
                        result.colUsed[localC]++;
                        result.usedCells++;

                        // Set boundary flags
                        if (r == r1) result.topBoundary[localC] = true;
                        if (r == r2) result.bottomBoundary[localC] = true;
                        if (c == c1) result.leftBoundary[localR] = true;
                        if (c == c2) result.rightBoundary[localR] = true;
                    }
                }
            }

            // Evaluate row constraints
            for (int r = r1; r <= r2; r++) {
                int used = countUsedInRow(r);
                if (used == rowClues[r]) {
                    score += 20;
                    result.exactRowMatches++;
                } else if (used < rowClues[r]) {
                    score += 5;
                } else {
                    score -= 25;
                }
            }

            // Evaluate column constraints
            for (int c = c1; c <= c2; c++) {
                int used = countUsedInCol(c);
                if (used == colClues[c]) {
                    score += 20;
                    result.exactColMatches++;
                } else if (used < colClues[c]) {
                    score += 5;
                } else {
                    score -= 25;
                }
            }

            // Evaluate path continuity
            int continuity = evaluateContinuityInQuadrant(r1, r2, c1, c2);
            score += continuity * 8;
            result.continuityScore = continuity;

            // Add closeness bonus to solver path
            int dist = manhattan(solverHead, playerHead);
            int closenessBonus = Math.max(0, 12 - dist) * 3;
            score += closenessBonus;

            result.score = score;
            return result;
        }

        /**
         * Combine four quadrants by synchronizing boundaries
         */
        private QuadrantResult combineQuadrants(int r1, int r2, int c1, int c2,
                                               QuadrantResult q1, QuadrantResult q2,
                                               QuadrantResult q3, QuadrantResult q4) {
            QuadrantResult parent = new QuadrantResult(r1, r2, c1, c2);
            int score = 0;
            int violations = 0;

            // Aggregate scores from quadrants
            if (q1 != null) score += q1.score;
            if (q2 != null) score += q2.score;
            if (q3 != null) score += q3.score;
            if (q4 != null) score += q4.score;

            // Check horizontal boundaries
            if (q1 != null && q2 != null) {
                int rows = Math.min(q1.rightBoundary.length, q2.leftBoundary.length);
                for (int i = 0; i < rows; i++) {
                    if (q1.rightBoundary[i] != q2.leftBoundary[i]) {
                        violations += 10;
                    } else if (q1.rightBoundary[i]) {
                        score += 8; // Good continuity
                    }
                }
            }

            if (q3 != null && q4 != null) {
                int rows = Math.min(q3.rightBoundary.length, q4.leftBoundary.length);
                for (int i = 0; i < rows; i++) {
                    if (q3.rightBoundary[i] != q4.leftBoundary[i]) {
                        violations += 10;
                    } else if (q3.rightBoundary[i]) {
                        score += 8;
                    }
                }
            }

            // Check vertical boundaries
            if (q1 != null && q3 != null) {
                int cols = Math.min(q1.bottomBoundary.length, q3.topBoundary.length);
                for (int i = 0; i < cols; i++) {
                    if (q1.bottomBoundary[i] != q3.topBoundary[i]) {
                        violations += 10;
                    } else if (q1.bottomBoundary[i]) {
                        score += 8;
                    }
                }
            }

            if (q2 != null && q4 != null) {
                int cols = Math.min(q2.bottomBoundary.length, q4.topBoundary.length);
                for (int i = 0; i < cols; i++) {
                    if (q2.bottomBoundary[i] != q4.topBoundary[i]) {
                        violations += 10;
                    } else if (q2.bottomBoundary[i]) {
                        score += 8;
                    }
                }
            }

            // Calculate aggregate row/col usage
            for (int r = r1; r <= r2; r++) {
                parent.rowUsed[r - r1] = countUsedInRow(r);
                parent.usedCells += parent.rowUsed[r - r1];
            }
            for (int c = c1; c <= c2; c++) {
                parent.colUsed[c - c1] = countUsedInCol(c);
            }

            // Set boundary flags for parent
            for (int c = c1; c <= c2; c++) {
                parent.topBoundary[c - c1] = grid[r1][c].isUsed();
                parent.bottomBoundary[c - c1] = grid[r2][c].isUsed();
            }
            for (int r = r1; r <= r2; r++) {
                parent.leftBoundary[r - r1] = grid[r][c1].isUsed();
                parent.rightBoundary[r - r1] = grid[r][c2].isUsed();
            }

            parent.boundaryViolations = violations;
            score -= violations;

            // Re-evaluate row/col constraints at parent level
            for (int r = r1; r <= r2; r++) {
                int used = countUsedInRow(r);
                if (used == rowClues[r]) {
                    score += 15;
                    parent.exactRowMatches++;
                } else if (used < rowClues[r]) {
                    score += 3;
                } else {
                    score -= 20;
                }
            }

            for (int c = c1; c <= c2; c++) {
                int used = countUsedInCol(c);
                if (used == colClues[c]) {
                    score += 15;
                    parent.exactColMatches++;
                } else if (used < colClues[c]) {
                    score += 3;
                } else {
                    score -= 20;
                }
            }

            // Bonus if heads are adjacent
            if (headsAdjacent()) score += 150;

            parent.score = score;
            return parent;
        }

        /**
         * Evaluate path continuity in a quadrant
         */
        private int evaluateContinuityInQuadrant(int r1, int r2, int c1, int c2) {
            int continuity = 0;
            
            for (int r = r1; r <= r2; r++) {
                for (int c = c1; c <= c2; c++) {
                    Cell cell = grid[r][c];
                    if (!cell.isUsed()) continue;

                    int neighbors = 0;
                    for (Cell n : getNeighbors(cell)) {
                        if (n.isUsed() && n.r >= r1 && n.r <= r2 && n.c >= c1 && n.c <= c2) {
                            neighbors++;
                        }
                    }

                    // Check if this cell has appropriate connections
                    if (neighbors == 1 && (cell == playerHead || cell == solverHead)) {
                        continuity += 2; // Good endpoint
                    } else if (neighbors == 2) {
                        continuity += 3; // Good mid-path
                    } else if (neighbors > 2) {
                        continuity -= 5; // Too many connections
                    } else if (neighbors == 0 && !cell.isUsed()) {
                        continuity -= 2; // Isolated used cell
                    }
                }
            }
            
            return continuity;
        }

        /**
         * Evaluate whole board score using D&C
         */
        private int evaluateWholeBoardScore() {
            QuadrantResult root = analyzeQuadrant(0, gridSize - 1, 0, gridSize - 1);
            return (root == null) ? Integer.MIN_VALUE / 2 : root.score;
        }

        /**
         * Choose best player move using D&C analysis
         */
        private Cell chooseBestPlayerMove() {
            List<Cell> validMoves = new ArrayList<>();
            for (Cell neighbor : getNeighbors(playerHead)) {
                if (isValidPlayerMove(neighbor.r, neighbor.c)) {
                    validMoves.add(neighbor);
                }
            }

            if (validMoves.isEmpty()) return null;

            Cell best = null;
            int bestScore = Integer.MIN_VALUE;

            for (Cell candidate : validMoves) {
                Cell oldHead = playerHead;
                candidate.fromPlayer = true;
                candidate.playerParent = playerHead;
                playerHead = candidate;

                int score = evaluateWholeBoardScore();
                score -= manhattan(playerHead, solverHead) * 2;
                if (headsAdjacent()) score += 200;

                undoPlayerMove(candidate, oldHead);

                if (score > bestScore) {
                    bestScore = score;
                    best = candidate;
                }
            }

            return best;
        }

        /**
         * Choose best solver move using D&C analysis
         * This is where the computer continues its path from B
         */
        private Cell chooseBestSolverMove() {
            List<Cell> validMoves = new ArrayList<>();
            for (Cell neighbor : getNeighbors(solverHead)) {
                if (isValidSolverMove(neighbor)) {
                    validMoves.add(neighbor);
                }
            }

            if (validMoves.isEmpty()) return null;

            Cell best = null;
            int bestScore = Integer.MIN_VALUE;

            for (Cell candidate : validMoves) {
                Cell oldHead = solverHead;
                candidate.fromSolver = true;
                candidate.solverParent = solverHead;
                solverHead = candidate;

                int score = evaluateWholeBoardScore();
                
                // Additional bonuses for solver path continuity
                score -= manhattan(solverHead, playerHead) * 2;
                score += evaluatePathContinuityFromSolver(candidate) * 10;
                
                if (headsAdjacent()) score += 200;

                undoSolverMove(candidate, oldHead);

                if (score > bestScore) {
                    bestScore = score;
                    best = candidate;
                }
            }

            return best;
        }

        /**
         * Evaluate continuity of solver's path
         */
        private int evaluatePathContinuityFromSolver(Cell candidate) {
            int continuity = 0;
            Cell current = candidate;
            int steps = 0;
            
            // Trace back along solver path
            while (current != null && steps < 5) {
                int neighbors = 0;
                for (Cell n : getNeighbors(current)) {
                    if (n.isOwnedBySolver()) neighbors++;
                }
                
                if (neighbors <= 2) continuity += 2;
                else continuity -= 3;
                
                current = current.solverParent;
                steps++;
            }
            
            return continuity;
        }

        /**
         * Computer's turn - continues its path from B
         */
        private void computerTurn() {
            if (gameOver) {
                computerThinking = false;
                return;
            }

            // Find valid moves for solver
            List<Cell> validMoves = new ArrayList<>();
            for (Cell neighbor : getNeighbors(solverHead)) {
                if (isValidSolverMove(neighbor)) validMoves.add(neighbor);
            }

            if (validMoves.isEmpty()) {
                computerThinking = false;
                gameOver = true;
                sidePanel.updateStatus("Game Over - Computer stuck");
                sidePanel.updateAnalysis("<html><font color='red'><b>No valid moves for computer!</b></font></html>");
                repaint();
                return;
            }

            // Use D&C to choose best move
            new Thread(() -> {
                Cell bestMove = chooseBestSolverMove();
                if (bestMove == null) bestMove = validMoves.get(0);

                final Cell move = bestMove;
                final int finalScore = evaluateWholeBoardScore();

                SwingUtilities.invokeLater(() -> {
                    computerThinking = false;
                    makeSolverMove(move);
                    
                    sidePanel.updateStatus("Computer moved to " + move);
                    sidePanel.updateAnalysis(
                        "<html><b>Quadrant D&C Analysis</b><br/>" +
                        "Computer continues from B → " + move + "<br/>" +
                        "Board Score: " + finalScore + "<br/>" +
                        "<b>Method:</b> 4-Way Boundary Sync</html>"
                    );
                    
                    repaint();
                    checkGameStatus();

                    if (!gameOver) {
                        isPlayerTurn = true;
                        sidePanel.updateStatus("Your Turn (A → ...)");
                    }
                });
            }).start();
        }

        /**
         * Auto-solver step
         */
        private void autoSolveStep() {
            if (!autoSolving || gameOver) {
                autoSolving = false;
                return;
            }

            // Player move
            Cell pMove = chooseBestPlayerMove();
            if (pMove == null) {
                autoSolving = false;
                gameOver = true;
                sidePanel.updateStatus("Auto-solver stopped");
                repaint();
                return;
            }

            makePlayerMove(pMove.r, pMove.c);
            repaint();
            checkGameStatus();
            
            if (gameOver) {
                autoSolving = false;
                return;
            }

            // Computer move
            Cell sMove = chooseBestSolverMove();
            if (sMove == null) {
                autoSolving = false;
                gameOver = true;
                sidePanel.updateStatus("Auto-solver stopped");
                repaint();
                return;
            }

            makeSolverMove(sMove);
            repaint();
            checkGameStatus();
            
            if (gameOver) {
                autoSolving = false;
                return;
            }

            sidePanel.updateAnalysis(
                "<html><b>Auto Solver</b><br/>" +
                "Player: " + pMove + "<br/>" +
                "Computer: " + sMove + "<br/>" +
                "Using Quadtree D&C</html>"
            );

            // Schedule next step
            new Timer().schedule(new TimerTask() {
                @Override
                public void run() {
                    SwingUtilities.invokeLater(() -> autoSolveStep());
                }
            }, 200);
        }

        /**
         * Start auto-solver
         */
        public void autoSolve() {
            if (gameOver || computerThinking || autoSolving) return;
            autoSolving = true;
            isPlayerTurn = false;
            clearHints();
            sidePanel.updateStatus("Auto-solving...");
            autoSolveStep();
        }

        /**
         * Show hint for player
         */
        public void showHint() {
            if (gameOver || computerThinking || autoSolving) return;
            clearHints();

            Cell best = chooseBestPlayerMove();

            if (best != null) {
                best.isHint = true;
                sidePanel.updateStatus("💡 Hint: Try " + best);
                repaint();
            } else {
                sidePanel.updateStatus("⚠ No valid moves!");
            }
        }

        /**
         * Clear hints
         */
        private void clearHints() {
            for (int r = 0; r < gridSize; r++) {
                for (int c = 0; c < gridSize; c++) {
                    grid[r][c].isHint = false;
                }
            }
        }

        /**
         * Check game status
         */
        private void checkGameStatus() {
            if (headsAdjacent()) {
                boolean valid = true;
                for (int i = 0; i < gridSize; i++) {
                    if (countUsedInRow(i) != rowClues[i] || countUsedInCol(i) != colClues[i]) {
                        valid = false;
                        break;
                    }
                }

                gameOver = true;
                autoSolving = false;

                if (valid) {
                    sidePanel.updateStatus("🎉 PUZZLE SOLVED! 🎉");
                    performAnalysis();
                } else {
                    sidePanel.updateStatus("❌ Paths met but clues unsatisfied");
                }
            }
        }

        /**
         * Perform final analysis
         */
        private void performAnalysis() {
            int totalSteps = totalCellsUsed();
            QuadrantResult finalResult = analyzeQuadrant(0, gridSize - 1, 0, gridSize - 1);

            String html = "<html><b>✓ Solution Verified!</b><br/>" +
                "<b>Algorithm:</b> Quadtree D&C<br/>" +
                "<b>Used Cells:</b> " + totalSteps + "<br/>" +
                "<b>Boundary Score:</b> " + finalResult.score + "<br/>" +
                "<b>Boundary Violations:</b> " + finalResult.boundaryViolations + "<br/>" +
                "<b>Continuity Score:</b> " + finalResult.continuityScore + "</html>";

            sidePanel.updateAnalysis(html);
        }

        /**
         * Paint the game board
         */
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int startX = 50, startY = 50;

            // Draw quadrant divider lines (visualizing D&C)
            g2.setColor(new Color(200, 200, 255, 100));
            g2.setStroke(new BasicStroke(2));
            int midX = startX + (gridSize / 2) * CELL_SIZE;
            int midY = startY + (gridSize / 2) * CELL_SIZE;
            g2.drawLine(midX, startY, midX, startY + gridSize * CELL_SIZE);
            g2.drawLine(startX, midY, startX + gridSize * CELL_SIZE, midY);

            // Draw grid and clues
            for (int r = 0; r < gridSize; r++) {
                // Row clues
                int rCur = countUsedInRow(r);
                Color rColor = (rCur == rowClues[r]) ? new Color(0, 150, 0)
                        : ((rCur > rowClues[r]) ? Color.RED : Color.BLACK);
                g2.setColor(rColor);
                g2.setFont(new Font("Monospaced", Font.BOLD, 14));
                g2.drawString(rCur + "/" + rowClues[r], startX - 45, startY + r * CELL_SIZE + 30);

                for (int c = 0; c < gridSize; c++) {
                    int x = startX + c * CELL_SIZE;
                    int y = startY + r * CELL_SIZE;

                    // Column clues
                    if (r == 0) {
                        int cCur = countUsedInCol(c);
                        Color cColor = (cCur == colClues[c]) ? new Color(0, 150, 0)
                                : ((cCur > colClues[c]) ? Color.RED : Color.BLACK);
                        g2.setColor(cColor);
                        g2.drawString(cCur + "/" + colClues[c], x + 12, startY - 10);
                    }

                    // Hint highlight
                    if (grid[r][c].isHint) {
                        g2.setColor(new Color(255, 255, 100, 200));
                        g2.fillRect(x + 1, y + 1, CELL_SIZE - 2, CELL_SIZE - 2);
                    }

                    // Grid lines
                    g2.setColor(Color.LIGHT_GRAY);
                    g2.setStroke(new BasicStroke(1));
                    g2.drawRect(x, y, CELL_SIZE, CELL_SIZE);

                    // Cross marks
                    if (grid[r][c].isCrossed) {
                        g2.setColor(new Color(255, 100, 100));
                        g2.setStroke(new BasicStroke(3));
                        g2.drawLine(x + 10, y + 10, x + CELL_SIZE - 10, y + CELL_SIZE - 10);
                        g2.drawLine(x + CELL_SIZE - 10, y + 10, x + 10, y + CELL_SIZE - 10);
                    }

                    // Start/End markers
                    if (grid[r][c] == startNode) {
                        g2.setColor(new Color(50, 100, 255));
                        g2.setFont(new Font("SansSerif", Font.BOLD, 20));
                        g2.drawString("A", x + 17, y + 33);
                    } else if (grid[r][c] == endNode) {
                        g2.setColor(new Color(255, 50, 50));
                        g2.setFont(new Font("SansSerif", Font.BOLD, 20));
                        g2.drawString("B", x + 17, y + 33);
                    }
                }
            }

            // Draw player path (green)
            g2.setColor(new Color(34, 139, 34));
            g2.setStroke(new BasicStroke(6, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));

            Cell tempP = playerHead;
            while (tempP != null && tempP.playerParent != null) {
                int x1 = 50 + tempP.c * CELL_SIZE + CELL_SIZE / 2;
                int y1 = 50 + tempP.r * CELL_SIZE + CELL_SIZE / 2;
                int x2 = 50 + tempP.playerParent.c * CELL_SIZE + CELL_SIZE / 2;
                int y2 = 50 + tempP.playerParent.r * CELL_SIZE + CELL_SIZE / 2;
                g2.drawLine(x1, y1, x2, y2);
                tempP = tempP.playerParent;
            }

            // Draw solver path (red)
            g2.setColor(new Color(180, 60, 60));
            g2.setStroke(new BasicStroke(6, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));

            Cell tempS = solverHead;
            while (tempS != null && tempS.solverParent != null) {
                int x1 = 50 + tempS.c * CELL_SIZE + CELL_SIZE / 2;
                int y1 = 50 + tempS.r * CELL_SIZE + CELL_SIZE / 2;
                int x2 = 50 + tempS.solverParent.c * CELL_SIZE + CELL_SIZE / 2;
                int y2 = 50 + tempS.solverParent.r * CELL_SIZE + CELL_SIZE / 2;
                g2.drawLine(x1, y1, x2, y2);
                tempS = tempS.solverParent;
            }

            // Draw connection if heads adjacent
            if (headsAdjacent()) {
                g2.setColor(new Color(120, 0, 180));
                g2.setStroke(new BasicStroke(4, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                int x1 = 50 + playerHead.c * CELL_SIZE + CELL_SIZE / 2;
                int y1 = 50 + playerHead.r * CELL_SIZE + CELL_SIZE / 2;
                int x2 = 50 + solverHead.c * CELL_SIZE + CELL_SIZE / 2;
                int y2 = 50 + solverHead.r * CELL_SIZE + CELL_SIZE / 2;
                g2.drawLine(x1, y1, x2, y2);
            }
        }
    }

    /**
     * Side panel with controls and information
     */
    class SidePanel extends JPanel {
        private JLabel statusLabel;
        private JLabel analysisLabel;

        public SidePanel(GamePanel gamePanel) {
            setPreferredSize(new Dimension(300, 640));
            setBackground(new Color(240, 240, 245));
            setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
            setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

            // Title
            JLabel title = new JLabel("TRACKS LAB");
            title.setFont(new Font("SansSerif", Font.BOLD, 28));
            title.setForeground(new Color(40, 80, 160));
            title.setAlignmentX(CENTER_ALIGNMENT);

            // Subtitle
            JLabel subtitle = new JLabel("Quadtree Divide & Conquer");
            subtitle.setFont(new Font("SansSerif", Font.ITALIC, 14));
            subtitle.setForeground(new Color(100, 100, 100));
            subtitle.setAlignmentX(CENTER_ALIGNMENT);

            // Instructions
            JLabel instr = new JLabel(
                "<html><center><b>HOW TO PLAY:</b><br/>" +
                "• Left Click: Add to your path from A<br/>" +
                "• Right Click: Mark X (block cell)<br/>" +
                "• Computer continues path from B<br/><br/>" +
                "<b>QUADTREE D&C:</b><br/>" +
                "• Divide board into 4 quadrants<br/>" +
                "• Analyze boundaries recursively<br/>" +
                "• Combine with boundary sync<br/>" +
                "• Computer chooses best continuation</center></html>"
            );
            instr.setAlignmentX(CENTER_ALIGNMENT);
            instr.setFont(new Font("SansSerif", Font.PLAIN, 11));

            // Status
            statusLabel = new JLabel("Your Turn (A → ...)");
            statusLabel.setForeground(new Color(0, 120, 180));
            statusLabel.setFont(new Font("SansSerif", Font.BOLD, 16));
            statusLabel.setAlignmentX(CENTER_ALIGNMENT);

            // Analysis area
            analysisLabel = new JLabel("<html><b>Ready</b></html>");
            analysisLabel.setAlignmentX(CENTER_ALIGNMENT);
            analysisLabel.setVerticalAlignment(SwingConstants.TOP);
            analysisLabel.setFont(new Font("Monospaced", Font.PLAIN, 11));
            analysisLabel.setBorder(BorderFactory.createLineBorder(new Color(200, 200, 200)));

            // Pattern selection
            JPanel patternPanel = new JPanel(new GridLayout(3, 2, 5, 5));
            patternPanel.setMaximumSize(new Dimension(250, 120));
            patternPanel.setBackground(new Color(240, 240, 245));

            for (int i = 1; i <= 5; i++) {
                final int pattern = i;
                JButton patBtn = new JButton("Pattern " + i);
                patBtn.setFont(new Font("SansSerif", Font.PLAIN, 11));
                patBtn.addActionListener(e -> gamePanel.loadPattern(pattern));
                patternPanel.add(patBtn);
            }

            JButton randomBtn = new JButton("Random");
            randomBtn.setFont(new Font("SansSerif", Font.PLAIN, 11));
            randomBtn.addActionListener(e -> gamePanel.newRandomPuzzle());
            patternPanel.add(randomBtn);

            // Action buttons
            JButton hintBtn = createStyledButton("💡 Show Hint", new Color(255, 200, 100));
            hintBtn.addActionListener(e -> gamePanel.showHint());

            JButton autoSolveBtn = createStyledButton("🤖 Auto Solve", new Color(100, 200, 255));
            autoSolveBtn.addActionListener(e -> gamePanel.autoSolve());

            JButton resetBtn = createStyledButton("🔄 Reset", new Color(200, 200, 200));
            resetBtn.addActionListener(e -> gamePanel.resetCurrentPuzzle());

            // Add all components
            add(title);
            add(Box.createRigidArea(new Dimension(0, 5)));
            add(subtitle);
            add(Box.createRigidArea(new Dimension(0, 15)));
            add(instr);
            add(Box.createRigidArea(new Dimension(0, 15)));
            add(statusLabel);
            add(Box.createRigidArea(new Dimension(0, 15)));

            JLabel patternLabel = new JLabel("Select Puzzle:");
            patternLabel.setAlignmentX(CENTER_ALIGNMENT);
            patternLabel.setFont(new Font("SansSerif", Font.BOLD, 14));
            add(patternLabel);
            add(Box.createRigidArea(new Dimension(0, 8)));
            add(patternPanel);

            add(Box.createRigidArea(new Dimension(0, 15)));
            add(hintBtn);
            add(Box.createRigidArea(new Dimension(0, 8)));
            add(autoSolveBtn);
            add(Box.createRigidArea(new Dimension(0, 8)));
            add(resetBtn);
            add(Box.createRigidArea(new Dimension(0, 15)));
            
            JLabel analysisTitle = new JLabel("Analysis:");
            analysisTitle.setAlignmentX(CENTER_ALIGNMENT);
            analysisTitle.setFont(new Font("SansSerif", Font.BOLD, 12));
            add(analysisTitle);
            add(Box.createRigidArea(new Dimension(0, 5)));
            add(analysisLabel);
        }

        /**
         * Create styled button
         */
        private JButton createStyledButton(String text, Color bgColor) {
            JButton button = new JButton(text);
            button.setAlignmentX(CENTER_ALIGNMENT);
            button.setMaximumSize(new Dimension(240, 40));
            button.setFont(new Font("SansSerif", Font.BOLD, 13));
            button.setBackground(bgColor);
            button.setOpaque(true);
            button.setBorderPainted(false);
            return button;
        }

        /**
         * Update status
         */
        public void updateStatus(String msg) {
            statusLabel.setText(msg);
        }

        /**
         * Update analysis
         */
        public void updateAnalysis(String html) {
            analysisLabel.setText(html);
        }
    }
}
