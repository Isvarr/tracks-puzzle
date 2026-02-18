import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;
import java.util.Timer;

/**
 * TracksGame - FINAL (50-50 Bidirectional Solver + Visible Forward/Backward Search)
 *
 * NEW:
 * - Forward explored cells = Light Blue
 * - Backward explored cells = Light Pink
 * - Shows solver is truly bidirectional
 */
public class TracksGame extends JFrame {

    private int gridSize = 8;
    private static final int CELL_SIZE = 50;

    private GamePanel gamePanel;
    private SidePanel sidePanel;
    private GameMode currentMode = GameMode.PATTERN1;

    enum GameMode { PATTERN1, PATTERN2, PATTERN3, PATTERN4, PATTERN5, RANDOM }

    public TracksGame() {
        setTitle("Arbiter Logic: Serpentine Solver (Bidirectional A* 50/50 + Visual Search)");
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
        try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); } catch (Exception ignored) {}
        SwingUtilities.invokeLater(() -> new TracksGame().setVisible(true));
    }

    static class Cell {
        int r, c;
        boolean hasTrack = false;
        boolean isCrossed = false;
        boolean isHint = false;
        Cell trackParent;

        public Cell(int r, int c) { this.r = r; this.c = c; }

        @Override
        public String toString() { return "(" + r + "," + c + ")"; }
    }

    class GamePanel extends JPanel {
        private Cell[][] grid;
        private int[] rowClues, colClues;
        private Cell startNode, endNode, currentHead;

        private boolean isPlayerTurn = true;
        private boolean gameOver = false;
        private boolean autoSolving = false;
        private boolean computerThinking = false;
        private List<Cell> solutionPath = null;

        // NEW: visual search highlight
        private BitSet forwardExplored = null;
        private BitSet backwardExplored = null;

        public GamePanel() {
            setBackground(Color.WHITE);
            initGame(GameMode.PATTERN1);

            addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent e) {
                    if (gameOver || !isPlayerTurn || autoSolving || computerThinking) return;

                    int startX = 50, startY = 50;
                    int c = (e.getX() - startX) / CELL_SIZE;
                    int r = (e.getY() - startY) / CELL_SIZE;

                    if (r >= 0 && c >= 0 && r < gridSize && c < gridSize) {

                        if (SwingUtilities.isRightMouseButton(e)) {
                            if (!grid[r][c].hasTrack && grid[r][c] != currentHead) {
                                grid[r][c].isCrossed = !grid[r][c].isCrossed;
                                repaint();
                            }
                            return;
                        }

                        if (SwingUtilities.isLeftMouseButton(e) && isValidMove(r, c)) {
                            clearHints();
                            makeMove(r, c);
                            repaint();
                            checkGameStatus();

                            if (!gameOver) {
                                computerThinking = true;
                                isPlayerTurn = false;
                                sidePanel.updateStatus("Computer thinking...");
                                new Timer().schedule(new TimerTask() {
                                    @Override public void run() { computerTurn(); }
                                }, 50);
                            }
                        }
                    }
                }
            });
        }

        public void initGame(GameMode mode) {
            currentMode = mode;
            switch(mode) {
                case PATTERN1: gridSize = 8; break;
                case PATTERN2: gridSize = 9; break;
                case PATTERN3: gridSize = 10; break;
                case PATTERN4: gridSize = 8; break;
                case PATTERN5: gridSize = 10; break;
                case RANDOM: gridSize = 8; break;
            }

            setPreferredSize(new Dimension(gridSize * CELL_SIZE + 100, gridSize * CELL_SIZE + 100));
            grid = new Cell[gridSize][gridSize];
            for (int r = 0; r < gridSize; r++)
                for (int c = 0; c < gridSize; c++)
                    grid[r][c] = new Cell(r, c);

            if (mode == GameMode.PATTERN1) {
                generateFixedPattern1();
            } else {
                long seed;
                switch(mode) {
                    case PATTERN2: seed = 12399L; break;
                    case PATTERN3: seed = 44455L; break;
                    case PATTERN4: seed = 99887L; break;
                    case PATTERN5: seed = 101010L; break;
                    default: seed = System.currentTimeMillis(); break;
                }
                generateHighComplexityPuzzle(seed);
            }

            currentHead = startNode;
            isPlayerTurn = true;
            gameOver = false;
            autoSolving = false;
            computerThinking = false;
            solutionPath = null;

            // clear search highlights
            forwardExplored = new BitSet(gridSize * gridSize);
            backwardExplored = new BitSet(gridSize * gridSize);

            TracksGame.this.pack();
            repaint();

            if (sidePanel != null) {
                String modeName = (mode == GameMode.RANDOM) ? "Random" : "Pattern " + mode.name().charAt(7);
                sidePanel.updateStatus("Player's Turn");
                sidePanel.updateAnalysis("<html><b>Puzzle:</b> " + modeName +
                        "<br/>Size: " + gridSize + "x" + gridSize +
                        "<br/>Length: " + sumArray(rowClues) + " cells</html>");
            }
        }

        private void generateFixedPattern1() {
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
            startNode.hasTrack = true;
            calculateCluesFromPath(path);
        }

        private void generateHighComplexityPuzzle(long seed) {
            Random rand = new Random(seed);
            int targetLength = (int)((gridSize * gridSize) * 0.50);

            while (true) {
                for(int r=0; r<gridSize; r++)
                    for(int c=0; c<gridSize; c++)
                        grid[r][c] = new Cell(r,c);

                int startR = rand.nextInt(gridSize - 2) + 1;
                int startC = rand.nextInt(gridSize - 2) + 1;
                Cell start = grid[startR][startC];
                boolean[][] visited = new boolean[gridSize][gridSize];
                visited[startR][startC] = true;

                List<Cell> bestPath = new ArrayList<>();
                List<Cell> currentPath = new ArrayList<>();
                currentPath.add(start);

                generatePathRecursive(start, visited, currentPath, bestPath, targetLength, rand);

                if (bestPath.size() >= targetLength) {
                    startNode = bestPath.get(0);
                    endNode = bestPath.get(bestPath.size()-1);
                    startNode.hasTrack = true;
                    calculateCluesFromPath(bestPath);
                    break;
                }
            }
        }

        private void generatePathRecursive(Cell current, boolean[][] visited, List<Cell> currentPath,
                                           List<Cell> bestPath, int target, Random rand) {
            if (currentPath.size() > bestPath.size()) {
                bestPath.clear();
                bestPath.addAll(currentPath);
            }
            if (bestPath.size() >= target) return;

            List<Cell> neighbors = getNeighbors(current);
            Collections.shuffle(neighbors, rand);

            for (Cell n : neighbors) {
                if (!visited[n.r][n.c] && countNeighbors(n, visited) <= 1) {
                    visited[n.r][n.c] = true;
                    currentPath.add(n);
                    generatePathRecursive(n, visited, currentPath, bestPath, target, rand);
                    if (bestPath.size() >= target) return;
                    currentPath.remove(currentPath.size() - 1);
                    visited[n.r][n.c] = false;
                }
            }
        }

        private void calculateCluesFromPath(List<Cell> path) {
            rowClues = new int[gridSize];
            colClues = new int[gridSize];
            for (Cell c : path) {
                rowClues[c.r]++;
                colClues[c.c]++;
            }
        }

        private int countNeighbors(Cell c, boolean[][] visited) {
            int count = 0;
            for (Cell n : getNeighbors(c)) if (visited[n.r][n.c]) count++;
            return count;
        }

        public void resetCurrentPuzzle() { initGame(currentMode); }
        public void newRandomPuzzle() { initGame(GameMode.RANDOM); }

        public void loadPattern(int patternNum) {
            switch(patternNum) {
                case 1: initGame(GameMode.PATTERN1); break;
                case 2: initGame(GameMode.PATTERN2); break;
                case 3: initGame(GameMode.PATTERN3); break;
                case 4: initGame(GameMode.PATTERN4); break;
                case 5: initGame(GameMode.PATTERN5); break;
            }
        }

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

        public boolean isValidMove(int r, int c) {
            if (r < 0 || c < 0 || r >= gridSize || c >= gridSize) return false;
            Cell target = grid[r][c];
            if (target.hasTrack || target.isCrossed) return false;
            if (Math.abs(r - currentHead.r) + Math.abs(c - currentHead.c) != 1) return false;
            if (countTrackInRow(r) >= rowClues[r]) return false;
            if (countTrackInCol(c) >= colClues[c]) return false;
            return true;
        }

        private void makeMove(int r, int c) {
            Cell nextCell = grid[r][c];
            nextCell.hasTrack = true;
            nextCell.trackParent = currentHead;
            currentHead = nextCell;
        }

        // --- MERGE SORT ---
        private void mergeSort(List<Cell> list, Comparator<Cell> comparator) {
            if (list.size() <= 1) return;
            int mid = list.size() / 2;
            List<Cell> left = new ArrayList<>(list.subList(0, mid));
            List<Cell> right = new ArrayList<>(list.subList(mid, list.size()));
            mergeSort(left, comparator);
            mergeSort(right, comparator);
            merge(list, left, right, comparator);
        }

        private void merge(List<Cell> result, List<Cell> left, List<Cell> right, Comparator<Cell> comp) {
            int i = 0, j = 0, k = 0;
            while (i < left.size() && j < right.size()) {
                if (comp.compare(left.get(i), right.get(j)) <= 0) result.set(k++, left.get(i++));
                else result.set(k++, right.get(j++));
            }
            while (i < left.size()) result.set(k++, left.get(i++));
            while (j < right.size()) result.set(k++, right.get(j++));
        }

        // ==========================================================
        // 50/50 BIDIRECTIONAL A* + VISUAL HIGHLIGHT
        // ==========================================================

        private int idx(int r, int c) { return r * gridSize + c; }

        private int manhattan(Cell a, Cell b) {
            return Math.abs(a.r - b.r) + Math.abs(a.c - b.c);
        }

        private String keyState(int cellId, int[] rowUsed, int[] colUsed) {
            StringBuilder sb = new StringBuilder(16 + gridSize * 6);
            sb.append(cellId).append('|');
            for (int i = 0; i < gridSize; i++) sb.append(rowUsed[i]).append(',');
            sb.append('|');
            for (int i = 0; i < gridSize; i++) sb.append(colUsed[i]).append(',');
            return sb.toString();
        }

        private boolean combineFeasible(int[] fRow, int[] fCol, int[] bRow, int[] bCol, Cell meet) {
            for (int i = 0; i < gridSize; i++) {
                int rTotal = fRow[i] + bRow[i] - ((meet.r == i) ? 1 : 0);
                if (rTotal != rowClues[i]) return false;

                int cTotal = fCol[i] + bCol[i] - ((meet.c == i) ? 1 : 0);
                if (cTotal != colClues[i]) return false;
            }
            return true;
        }

        private boolean disjointExceptMeet(BitSet fVis, BitSet bVis, Cell meet) {
            BitSet tmp = (BitSet) fVis.clone();
            tmp.and(bVis);
            int meetIndex = idx(meet.r, meet.c);
            if (tmp.isEmpty()) return true;
            return tmp.cardinality() == 1 && tmp.get(meetIndex);
        }

        class FState {
            Cell pos;
            int g, f;
            int[] rowUsed, colUsed;
            BitSet visited;
            FState parent;

            FState(Cell pos, int g, int f, int[] rowUsed, int[] colUsed, BitSet visited, FState parent) {
                this.pos = pos; this.g = g; this.f = f;
                this.rowUsed = rowUsed; this.colUsed = colUsed;
                this.visited = visited; this.parent = parent;
            }
        }

        class BState {
            Cell pos;
            int g, f;
            int[] rowUsed, colUsed;
            BitSet visited;
            BState parent;

            BState(Cell pos, int g, int f, int[] rowUsed, int[] colUsed, BitSet visited, BState parent) {
                this.pos = pos; this.g = g; this.f = f;
                this.rowUsed = rowUsed; this.colUsed = colUsed;
                this.visited = visited; this.parent = parent;
            }
        }

        private List<Cell> reconstructBidirectional(FState fMeet, BState bMeet, Cell forwardStart) {
            ArrayList<Cell> left = new ArrayList<>();
            FState ft = fMeet;
            while (ft != null) { left.add(ft.pos); ft = ft.parent; }
            Collections.reverse(left);

            ArrayList<Cell> right = new ArrayList<>();
            BState bt = bMeet.parent;
            while (bt != null) { right.add(bt.pos); bt = bt.parent; }

            left.addAll(right);

            if (!left.isEmpty() && left.get(0) == forwardStart) left.remove(0);
            return left;
        }

        private long getSmartTimeLimitMs() {
            if (gridSize <= 8) return 7000;
            if (gridSize == 9) return 9000;
            return 18000;
        }

        private void markExploredForward(Cell c) {
            if (forwardExplored == null) return;
            forwardExplored.set(idx(c.r, c.c));
        }

        private void markExploredBackward(Cell c) {
            if (backwardExplored == null) return;
            backwardExplored.set(idx(c.r, c.c));
        }

        private List<Cell> runBidirectionalAStar(Cell forwardStart) {

            long startTime = System.currentTimeMillis();
            long timeLimitMs = getSmartTimeLimitMs();

            BitSet fixed = new BitSet(gridSize * gridSize);
            int[] fRow0 = new int[gridSize];
            int[] fCol0 = new int[gridSize];

            Cell t = forwardStart;
            while (t != null) {
                int id = idx(t.r, t.c);
                if (!fixed.get(id)) {
                    fixed.set(id);
                    fRow0[t.r]++;
                    fCol0[t.c]++;
                }
                t = t.trackParent;
            }

            int totalNeeded = sumArray(rowClues);
            int alreadyUsed = 0;
            for (int v : fRow0) alreadyUsed += v;
            int remainingBudget = totalNeeded - alreadyUsed;
            if (remainingBudget < 0) return null;

            FState fStart = new FState(
                    forwardStart,
                    0,
                    manhattan(forwardStart, endNode),
                    fRow0,
                    fCol0,
                    (BitSet) fixed.clone(),
                    null
            );

            BitSet bVis0 = new BitSet(gridSize * gridSize);
            bVis0.set(idx(endNode.r, endNode.c));
            int[] bRow0 = new int[gridSize];
            int[] bCol0 = new int[gridSize];
            bRow0[endNode.r] = 1;
            bCol0[endNode.c] = 1;

            BState bStart = new BState(
                    endNode,
                    0,
                    manhattan(endNode, forwardStart),
                    bRow0,
                    bCol0,
                    bVis0,
                    null
            );

            PriorityQueue<FState> fOpen = new PriorityQueue<>(Comparator.comparingInt(s -> s.f));
            PriorityQueue<BState> bOpen = new PriorityQueue<>(Comparator.comparingInt(s -> s.f));
            fOpen.add(fStart);
            bOpen.add(bStart);

            HashSet<String> fClosed = new HashSet<>();
            HashSet<String> bClosed = new HashSet<>();

            ArrayList<ArrayList<FState>> fAt = new ArrayList<>();
            ArrayList<ArrayList<BState>> bAt = new ArrayList<>();
            for (int i = 0; i < gridSize * gridSize; i++) {
                fAt.add(new ArrayList<>());
                bAt.add(new ArrayList<>());
            }
            fAt.get(idx(fStart.pos.r, fStart.pos.c)).add(fStart);
            bAt.get(idx(bStart.pos.r, bStart.pos.c)).add(bStart);

            int maxIterations = 600000;
            int it = 0;

            boolean expandForward = false; // so first flip becomes forward

            while (!fOpen.isEmpty() && !bOpen.isEmpty()) {

                if (System.currentTimeMillis() - startTime > timeLimitMs) return null;
                if (++it > maxIterations) break;

                expandForward = !expandForward;

                if (expandForward) {
                    FState cur = fOpen.poll();
                    if (cur == null) continue;
                    if (cur.g > remainingBudget) continue;

                    // VISUAL
                    markExploredForward(cur.pos);
                    if (it % 1200 == 0) repaint();

                    String curKey = keyState(idx(cur.pos.r, cur.pos.c), cur.rowUsed, cur.colUsed);
                    if (!fClosed.add(curKey)) continue;

                    int meetId = idx(cur.pos.r, cur.pos.c);
                    for (BState other : bAt.get(meetId)) {
                        if (cur.g + other.g > remainingBudget) continue;
                        if (!disjointExceptMeet(cur.visited, other.visited, cur.pos)) continue;
                        if (!combineFeasible(cur.rowUsed, cur.colUsed, other.rowUsed, other.colUsed, cur.pos)) continue;
                        return reconstructBidirectional(cur, other, forwardStart);
                    }

                    for (Cell nxt : getNeighbors(cur.pos)) {
                        int nid = idx(nxt.r, nxt.c);
                        if (grid[nxt.r][nxt.c].isCrossed) continue;
                        if (cur.visited.get(nid)) continue;

                        int[] nr = Arrays.copyOf(cur.rowUsed, gridSize);
                        int[] nc = Arrays.copyOf(cur.colUsed, gridSize);
                        nr[nxt.r]++; nc[nxt.c]++;
                        if (nr[nxt.r] > rowClues[nxt.r]) continue;
                        if (nc[nxt.c] > colClues[nxt.c]) continue;

                        BitSet nVis = (BitSet) cur.visited.clone();
                        nVis.set(nid);

                        int ng = cur.g + 1;
                        int nf = ng + manhattan(nxt, endNode);

                        FState ns = new FState(nxt, ng, nf, nr, nc, nVis, cur);

                        String nk = keyState(nid, nr, nc);
                        if (fClosed.contains(nk)) continue;

                        fAt.get(nid).add(ns);
                        fOpen.add(ns);
                    }

                } else {
                    BState cur = bOpen.poll();
                    if (cur == null) continue;
                    if (cur.g > remainingBudget) continue;

                    // VISUAL
                    markExploredBackward(cur.pos);
                    if (it % 1200 == 0) repaint();

                    String curKey = keyState(idx(cur.pos.r, cur.pos.c), cur.rowUsed, cur.colUsed);
                    if (!bClosed.add(curKey)) continue;

                    int meetId = idx(cur.pos.r, cur.pos.c);
                    for (FState other : fAt.get(meetId)) {
                        if (cur.g + other.g > remainingBudget) continue;
                        if (!disjointExceptMeet(other.visited, cur.visited, cur.pos)) continue;
                        if (!combineFeasible(other.rowUsed, other.colUsed, cur.rowUsed, cur.colUsed, cur.pos)) continue;
                        return reconstructBidirectional(other, cur, forwardStart);
                    }

                    for (Cell nxt : getNeighbors(cur.pos)) {
                        int nid = idx(nxt.r, nxt.c);
                        if (grid[nxt.r][nxt.c].isCrossed) continue;
                        if (fixed.get(nid)) continue;
                        if (cur.visited.get(nid)) continue;

                        int[] nr = Arrays.copyOf(cur.rowUsed, gridSize);
                        int[] nc = Arrays.copyOf(cur.colUsed, gridSize);
                        nr[nxt.r]++; nc[nxt.c]++;
                        if (nr[nxt.r] > rowClues[nxt.r]) continue;
                        if (nc[nxt.c] > colClues[nxt.c]) continue;

                        BitSet nVis = (BitSet) cur.visited.clone();
                        nVis.set(nid);

                        int ng = cur.g + 1;
                        int nf = ng + manhattan(nxt, forwardStart);

                        BState ns = new BState(nxt, ng, nf, nr, nc, nVis, cur);

                        String nk = keyState(nid, nr, nc);
                        if (bClosed.contains(nk)) continue;

                        bAt.get(nid).add(ns);
                        bOpen.add(ns);
                    }
                }
            }

            return null;
        }

        private List<Cell> findSolution() {
            return runBidirectionalAStar(currentHead);
        }

        private void gracefulSolverExit(String title, String reason) {
            autoSolving = false;
            computerThinking = false;
            isPlayerTurn = false;
            gameOver = true;

            sidePanel.updateStatus(title);
            sidePanel.updateAnalysis(
                    "<html><b>" + reason + "</b><br/><br/>Try another pattern or click <b>Reset</b>.</html>"
            );
            repaint();
        }

        // --- COMPUTER TURN ---
        private void computerTurn() {
            if (gameOver) { computerThinking = false; return; }

            List<Cell> validMoves = new ArrayList<>();
            for (Cell neighbor : getNeighbors(currentHead)) {
                if (isValidMove(neighbor.r, neighbor.c)) validMoves.add(neighbor);
            }

            if (validMoves.isEmpty()) {
                computerThinking = false;
                gameOver = true;
                sidePanel.updateStatus("Game Over - Draw");
                sidePanel.updateAnalysis("<html><font color='orange'><b>No Moves!</b></font></html>");
                repaint();
                return;
            }

            mergeSort(validMoves, (a, b) -> {
                int distA = Math.abs(a.r - endNode.r) + Math.abs(a.c - endNode.c);
                int distB = Math.abs(b.r - endNode.r) + Math.abs(b.c - endNode.c);
                return Integer.compare(distA, distB);
            });

            new Thread(() -> {
                Cell bestMove = null;
                for (Cell candidate : validMoves) {
                    Cell originalParent = candidate.trackParent;
                    candidate.trackParent = currentHead;

                    if (runBidirectionalAStar(candidate) != null) {
                        bestMove = candidate;
                        candidate.trackParent = originalParent;
                        break;
                    }
                    candidate.trackParent = originalParent;
                }

                final Cell move = bestMove;
                SwingUtilities.invokeLater(() -> {
                    computerThinking = false;
                    if (move != null) {
                        makeMove(move.r, move.c);
                        sidePanel.updateStatus("Computer → " + move);
                    } else {
                        gracefulSolverExit("⚠ Solver Stuck",
                                "Computer found no safe move that guarantees solvability.");
                        return;
                    }
                    repaint();
                    checkGameStatus();
                    if (!gameOver) { isPlayerTurn = true; sidePanel.updateStatus("Your Turn"); }
                });
            }).start();
        }

        public void showHint() {
            if (gameOver || autoSolving || computerThinking) return;
            clearHints();

            List<Cell> validMoves = new ArrayList<>();
            for (Cell neighbor : getNeighbors(currentHead)) {
                if (isValidMove(neighbor.r, neighbor.c)) validMoves.add(neighbor);
            }

            if (validMoves.isEmpty()) { sidePanel.updateStatus("⚠ No valid moves!"); return; }

            mergeSort(validMoves, (a, b) -> {
                int distA = Math.abs(a.r - endNode.r) + Math.abs(a.c - endNode.c);
                int distB = Math.abs(b.r - endNode.r) + Math.abs(b.c - endNode.c);
                return Integer.compare(distA, distB);
            });

            new Thread(() -> {
                Cell hintMove = null;
                for (Cell candidate : validMoves) {
                    Cell originalParent = candidate.trackParent;
                    candidate.trackParent = currentHead;

                    if (runBidirectionalAStar(candidate) != null) {
                        hintMove = candidate;
                        candidate.trackParent = originalParent;
                        break;
                    }
                    candidate.trackParent = originalParent;
                }

                final Cell finalHint = hintMove;
                SwingUtilities.invokeLater(() -> {
                    if (finalHint != null) {
                        grid[finalHint.r][finalHint.c].isHint = true;
                        sidePanel.updateStatus("💡 Hint: " + finalHint);
                    } else {
                        sidePanel.updateStatus("⚠ No solution found!");
                    }
                    repaint();
                });
            }).start();
        }

        private void clearHints() {
            for (int r = 0; r < gridSize; r++)
                for (int c = 0; c < gridSize; c++)
                    grid[r][c].isHint = false;
        }

        public void autoSolve() {
            if (gameOver || autoSolving || computerThinking) return;
            autoSolving = true;
            clearHints();

            // reset highlights every solve
            forwardExplored.clear();
            backwardExplored.clear();

            sidePanel.updateStatus("Auto-solving...");

            new Thread(() -> {
                List<Cell> solution = findSolution();
                SwingUtilities.invokeLater(() -> {
                    if (solution != null && !solution.isEmpty()) {
                        solutionPath = solution;
                        animateSolution(0);
                    } else {
                        gracefulSolverExit("⚠ Solver Timed Out",
                                "Solver could not finish within the time limit (puzzle may still be solvable).");
                    }
                });
            }).start();
        }

        private void animateSolution(int index) {
            if (solutionPath == null || index >= solutionPath.size()) {
                autoSolving = false;
                checkGameStatus();
                return;
            }
            Cell nextCell = solutionPath.get(index);
            if (isValidMove(nextCell.r, nextCell.c)) {
                makeMove(nextCell.r, nextCell.c);
                repaint();
                new Timer().schedule(new TimerTask() {
                    @Override public void run() {
                        SwingUtilities.invokeLater(() -> animateSolution(index + 1));
                    }
                }, 100);
            } else {
                autoSolving = false;
                gracefulSolverExit("⚠ Auto-solve stopped",
                        "Auto-solve stopped due to invalid step (stale solution).");
            }
        }

        private int countTrackInRow(int r) {
            int count = 0;
            for (int c = 0; c < gridSize; c++) if (grid[r][c].hasTrack) count++;
            return count;
        }

        private int countTrackInCol(int c) {
            int count = 0;
            for (int r = 0; r < gridSize; r++) if (grid[r][c].hasTrack) count++;
            return count;
        }

        private int sumArray(int[] arr) {
            int sum = 0;
            for (int v : arr) sum += v;
            return sum;
        }

        private void checkGameStatus() {
            if (currentHead == endNode) {
                gameOver = true;
                boolean valid = true;
                for (int i = 0; i < gridSize; i++) {
                    if (countTrackInRow(i) != rowClues[i] || countTrackInCol(i) != colClues[i]) {
                        valid = false; break;
                    }
                }
                if (valid) {
                    sidePanel.updateStatus("🎉 PUZZLE SOLVED! 🎉");
                    performAnalysis();
                } else {
                    sidePanel.updateStatus("❌ Constraints Failed");
                }
            }
        }

        private void performAnalysis() {
            int totalSteps = 0;
            Cell temp = currentHead;
            while (temp != null) { totalSteps++; temp = temp.trackParent; }
            String html = "<html><font color='black'>" +
                        "<b>✓ Solution Verified!</b><br/>" +
                        "<b>Algorithm:</b> Bidirectional A* (50/50) + MergeSort<br/>" +
                        "<b>Steps:</b> " + totalSteps + "<br/>" +
                        "<b>Coverage:</b> " + ((totalSteps * 100) / (gridSize * gridSize)) + "%" +
                        "</font></html>"
;
            sidePanel.updateAnalysis(html);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int startX = 50, startY = 50;

            for (int r = 0; r < gridSize; r++) {
                g2.setFont(new Font("SansSerif", Font.BOLD, 14));
                int rCur = countTrackInRow(r);
                Color rColor = (rCur == rowClues[r]) ? new Color(0, 150, 0) :
                        ((rCur > rowClues[r]) ? Color.RED : Color.BLACK);
                g2.setColor(rColor);
                g2.drawString(rCur + "/" + rowClues[r], startX - 45, startY + r * CELL_SIZE + 30);

                for (int c = 0; c < gridSize; c++) {
                    int x = startX + c * CELL_SIZE;
                    int y = startY + r * CELL_SIZE;

                    if (r == 0) {
                        int cCur = countTrackInCol(c);
                        Color cColor = (cCur == colClues[c]) ? new Color(0, 150, 0) :
                                ((cCur > colClues[c]) ? Color.RED : Color.BLACK);
                        g2.setColor(cColor);
                        g2.drawString(cCur + "/" + colClues[c], x + 12, startY - 10);
                    }

                    // NEW: show explored cells
                    int id = idx(r, c);
                    if (forwardExplored != null && forwardExplored.get(id)) {
                        g2.setColor(new Color(150, 200, 255, 120)); // light blue
                        g2.fillRect(x + 1, y + 1, CELL_SIZE - 2, CELL_SIZE - 2);
                    }
                    if (backwardExplored != null && backwardExplored.get(id)) {
                        g2.setColor(new Color(255, 170, 190, 120)); // light pink
                        g2.fillRect(x + 1, y + 1, CELL_SIZE - 2, CELL_SIZE - 2);
                    }

                    if (grid[r][c].isHint) {
                        g2.setColor(new Color(255, 255, 100, 200));
                        g2.fillRect(x + 1, y + 1, CELL_SIZE - 2, CELL_SIZE - 2);
                    }

                    g2.setColor(Color.LIGHT_GRAY);
                    g2.setStroke(new BasicStroke(1));
                    g2.drawRect(x, y, CELL_SIZE, CELL_SIZE);

                    if (grid[r][c].isCrossed) {
                        g2.setColor(new Color(255, 100, 100));
                        g2.setStroke(new BasicStroke(3));
                        g2.drawLine(x+10, y+10, x+CELL_SIZE-10, y+CELL_SIZE-10);
                        g2.drawLine(x+CELL_SIZE-10, y+10, x+10, y+CELL_SIZE-10);
                    }

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

            g2.setColor(new Color(34, 139, 34));
            g2.setStroke(new BasicStroke(6, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));

            Cell temp = currentHead;
            while (temp != null && temp.trackParent != null) {
                int x1 = startX + temp.c * CELL_SIZE + CELL_SIZE/2;
                int y1 = startY + temp.r * CELL_SIZE + CELL_SIZE/2;
                int x2 = startX + temp.trackParent.c * CELL_SIZE + CELL_SIZE/2;
                int y2 = startY + temp.trackParent.r * CELL_SIZE + CELL_SIZE/2;
                g2.drawLine(x1, y1, x2, y2);
                temp = temp.trackParent;
            }
        }
    }

    class SidePanel extends JPanel {
        private JLabel statusLabel;
        private JLabel analysisLabel;
        private GamePanel gamePanel;

        public SidePanel(GamePanel gamePanel) {
            this.gamePanel = gamePanel;
            setPreferredSize(new Dimension(260, 600));
            setBackground(new Color(240, 240, 245));
            setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
            setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

            JLabel title = new JLabel("Tracks Lab");
            title.setFont(new Font("SansSerif", Font.BOLD, 26));
            title.setForeground(new Color(40, 80, 160));
            title.setAlignmentX(CENTER_ALIGNMENT);

            JLabel instr = new JLabel("<html><center><font color='black'><b>Controls:</b><br/>Left Click: Add Track<br/>Right Click: Mark X<br/><br/><b>Goal:</b> Connect A→B<br/>Match All Clues</font></center></html>");
            instr.setAlignmentX(CENTER_ALIGNMENT);
            instr.setFont(new Font("SansSerif", Font.PLAIN, 12));

            statusLabel = new JLabel("Player's Turn");
            statusLabel.setForeground(new Color(0, 120, 180));
            statusLabel.setFont(new Font("SansSerif", Font.BOLD, 16));
            statusLabel.setAlignmentX(CENTER_ALIGNMENT);

            analysisLabel = new JLabel("");
            analysisLabel.setAlignmentX(CENTER_ALIGNMENT);
            analysisLabel.setVerticalAlignment(SwingConstants.TOP);
            analysisLabel.setFont(new Font("SansSerif", Font.PLAIN, 11));

            JPanel patternPanel = new JPanel();
            patternPanel.setLayout(new GridLayout(3, 2, 5, 5));
            patternPanel.setMaximumSize(new Dimension(230, 120));
            patternPanel.setBackground(new Color(240, 240, 245));

            for (int i = 1; i <= 5; i++) {
                final int pattern = i;
                JButton patBtn = new JButton("Pattern " + i);
                patBtn.setFont(new Font("SansSerif", Font.PLAIN, 11));
                patBtn.addActionListener(e -> gamePanel.loadPattern(pattern));
                patBtn.setForeground(Color.BLACK);
                patternPanel.add(patBtn);
            }

            JButton randomBtn = new JButton("Random");
            randomBtn.setFont(new Font("SansSerif", Font.PLAIN, 11));
            randomBtn.addActionListener(e -> gamePanel.newRandomPuzzle());
            randomBtn.setForeground(Color.BLACK);
            patternPanel.add(randomBtn);

            JButton hintBtn = new JButton("💡 Show Hint");
            hintBtn.setAlignmentX(CENTER_ALIGNMENT);
            hintBtn.setMaximumSize(new Dimension(210, 38)); 
            hintBtn.setFont(new Font("SansSerif", Font.BOLD, 13));
            hintBtn.setForeground(Color.BLACK);
            hintBtn.addActionListener(e -> gamePanel.showHint());

            JButton autoSolveBtn = new JButton("🤖 Auto-Solve");
            autoSolveBtn.setAlignmentX(CENTER_ALIGNMENT);
            autoSolveBtn.setMaximumSize(new Dimension(210, 38));
            autoSolveBtn.setFont(new Font("SansSerif", Font.BOLD, 13));
            autoSolveBtn.setForeground(Color.BLACK);
            autoSolveBtn.addActionListener(e -> gamePanel.autoSolve());

            JButton resetBtn = new JButton("🔄 Reset");
            resetBtn.setAlignmentX(CENTER_ALIGNMENT);
            resetBtn.setMaximumSize(new Dimension(210, 38));
            resetBtn.setFont(new Font("SansSerif", Font.BOLD, 13));
            resetBtn.setForeground(Color.BLACK);
            resetBtn.addActionListener(e -> gamePanel.resetCurrentPuzzle());

            add(title);
            add(Box.createRigidArea(new Dimension(0, 12)));
            add(instr);
            add(Box.createRigidArea(new Dimension(0, 20)));
            add(statusLabel);
            add(Box.createRigidArea(new Dimension(0, 18)));

            JLabel patternLabel = new JLabel("Select Puzzle:");
            patternLabel.setAlignmentX(CENTER_ALIGNMENT);
            patternLabel.setFont(new Font("SansSerif", Font.BOLD, 14));
            patternLabel.setForeground(Color.BLACK);
            add(patternLabel);
            add(Box.createRigidArea(new Dimension(0, 8)));
            add(patternPanel);

            add(Box.createRigidArea(new Dimension(0, 18)));
            add(hintBtn);
            add(Box.createRigidArea(new Dimension(0, 10)));
            add(autoSolveBtn);
            add(Box.createRigidArea(new Dimension(0, 10)));
            add(resetBtn);
            add(Box.createRigidArea(new Dimension(0, 25)));
            add(analysisLabel);
        }

        public void updateStatus(String msg) { statusLabel.setText(msg); }
        public void updateAnalysis(String html) { analysisLabel.setText(html); }
    }
}
