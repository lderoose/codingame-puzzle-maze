#include <set>
#include <queue>
#include <chrono>
#include <string>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <stdexcept>

using namespace std;


class Timer
{
    private:
        chrono::high_resolution_clock::time_point start_time;
        chrono::high_resolution_clock::time_point stop_time;
    public:
        void start() {start_time = chrono::high_resolution_clock::now();}
        void stop() {stop_time = chrono::high_resolution_clock::now();}
        double measure() {return chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start_time).count();}
};

class Node
{
    public:
        int row, col;
        string value = "?";
        pair<int, int> pos_neighbors[4];  // order: U, R, D, L
        bool visited = false;
        bool neighbors_filled = false;
        Node* parent = nullptr;
        Node(){}
        Node(int row, int col) : row(row), col(col) {}
        ~Node(){}
        Node(const Node& n) : row(n.row), col(n.col), value(n.value), visited(n.visited), neighbors_filled(n.neighbors_filled) 
        {
            for (int i=0; i<4; i++)
            {
                pos_neighbors[i].first = n.pos_neighbors[i].first;
                pos_neighbors[i].second = n.pos_neighbors[i].second;
            }
        }
        inline bool operator==(const Node& n) {return row==n.row && col==n.col;}
        pair<int, int> to_pair() {return make_pair(row, col);}
};

class Maze
{
    public:
        int n_rows;
        int n_cols;
        Node **grid;
        pair<int, int> pos_command = make_pair(-1, -1);
        bool command_found = false;
        bool command_reachable = false;
        bool alarm_enable = false;
        set<pair<int, int>> walkable;
        Maze(){};
        Maze(int n_rows, int n_cols) : n_rows(n_rows), n_cols(n_cols)
        {
            grid = new Node*[n_rows];
            for (int i=0; i<n_rows; i++) 
            {
                grid[i] = new Node[n_cols];
            }
        }
        ~Maze()
        {
            for (int i=0; i<n_rows; i++) {delete[] grid[i];}
            delete[] grid;
        }
        inline Node* operator[](int row){return grid[row];}
        inline string getValue(int row, int col) const {return grid[row][col].value;}
        inline void findNeighbors(int row, int col)
        {
            bool fill[4] = {false, false, false, false};
            Node n_up = grid[max(0, row-1)][col];
            Node n_right = grid[row][min(n_cols, col+1)];
            Node n_down = grid[min(n_rows, row+1)][col];
            Node n_left = grid[row][max(0, col-1)];
            if (n_up.value != "?") 
            {
                grid[row][col].pos_neighbors[0].first = n_up.row;
                grid[row][col].pos_neighbors[0].second = n_up.col;
                fill[0] = true;
            }
            if (n_right.value != "?") 
            {
                grid[row][col].pos_neighbors[1].first = n_right.row;
                grid[row][col].pos_neighbors[1].second = n_right.col;
                fill[1] = true;
            }
            if (n_down.value != "?") 
            {
                grid[row][col].pos_neighbors[2].first = n_down.row;
                grid[row][col].pos_neighbors[2].second = n_down.col;
                fill[2] = true;
            }
            if (n_left.value != "?") 
            {
                grid[row][col].pos_neighbors[3].first = n_left.row;
                grid[row][col].pos_neighbors[3].second = n_left.col;
                fill[3] = true;
            }

            // at border no neighbors
            if (row == 0) {fill[0] = true;}
            if (col == (n_cols-1)) {fill[1] = true;}
            if (row == (n_rows-1)) {fill[2] = true;}
            if (col == 0) {fill[3] = true;}

            // check if neighbors are totally filled
            grid[row][col].neighbors_filled = true;  // init to true
            for (int i=0; i<4; i++)
            {
                if (!fill[i]) 
                {
                    grid[row][col].neighbors_filled = false;
                    break;
                }
            }
        }
        void updateData(int iteration)
        {
            string row; // C of the characters in '#.TC?' (i.e. one line of the ASCII maze).
            for (int i=0; i<n_rows; i++) 
            {
                // read
                cin >> row; cin.ignore();

                // update maze values
                for (int j=0; j<n_cols; j++)
                {
                    char value = row[j];
                    grid[i][j].value = value;
                    grid[i][j].visited = false; // reset visited

                    // update node value row col one time
                    if (iteration == 0)
                    {
                        grid[i][j].row = i;
                        grid[i][j].col = j;
                    }

                    // update node if not already done after iteration 0
                    else
                    {
                        if (grid[i][j].neighbors_filled == false) {findNeighbors(i, j);}
                    }

                    if (value!='#' && value !='?') {walkable.insert(make_pair(i, j));}

                    // update if  command is found
                    if (value == 'C') 
                    {
                        pos_command.first = i;
                        pos_command.second = j;
                        command_found = true;
                    }
                }
            }
        }  
};


class Agent
{
    public: 
        int row = -1;
        int col = -1;
        int row_start = row;
        int col_start = col;
        set<pair<int, int>> visited;
        Agent(){}
        Agent(int row, int col): row(row), col(col), row_start(row), col_start(col){};
        ~Agent(){}
        inline pair<int, int> getPosition() const{return make_pair(row, col);}
        void initFirstPosition(int iteration, int new_row, int new_col)
        {
            if (iteration == 0)
            {
                row = new_row;
                col = new_col;
                row_start = row;
                col_start = col;
            }
        }
        inline void move(string &action)
        // the validity of the movement is supposed to be realized upstream.
        {
            if (action == "UP") {this->row -= 1;}
            else if (action == "DOWN") {this->row += 1;}
            else if (action == "LEFT") {this->col -= 1;}
            else if (action == "RIGHT") {this->col += 1;}
            else 
            {
                string fname = __FUNCTION__;
                throw invalid_argument(fname + ": action must belongs to UP, DOWN, LEFT, RIGHT, get " + action + " instead.");
            }
            visited.insert(make_pair(row, col));
        }
};


bool bfs(Node* start, Node* end, Maze* maze, vector<Node*>& nodesPath)
{
    queue<Node*> q;
    start->visited = true;
    q.push(start);

    while (!q.empty()) 
    {
        Node* current = q.front();
        q.pop();

        if (*current == *end) 
        {
            while (current) 
            {
                if (current->visited == false) {break;}
                current->visited = false;
                nodesPath.push_back(current);
                Node* old_current = current;
                current = current->parent;
            }
            return true;
        }    

        for (pair<int, int> pos_neighbor : (*current).pos_neighbors)
        {
            Node* neighbor = maze[pos_neighbor.first][pos_neighbor.second];
            if (!neighbor->visited && neighbor) 
            {
                if (neighbor->value!="#" && neighbor->value!="?")
                {
                    neighbor->visited = true;
                    neighbor->parent = current;
                    q.push(neighbor);
                }
            }
        }
    }
    return false;
}

bool isReachable(Node* start, Node* end, Maze* maze)
{
    bool reachable = false;
    vector<Node*> nodesPath;
    if (start->value == "C") {reachable = bfs(start, end, maze, nodesPath);}
    return reachable;
}


string findMove(Agent* agent, Node* node)
{
    int node_row = node->row;
    int node_col = node->col;
    int agent_row = agent->getPosition().first;
    int agent_col = agent->getPosition().second;

    string move = "";
    if (node_row < agent_row) {move = "UP";}
    if (node_row > agent_row) {move = "DOWN";}
    if (node_col < agent_col) {move = "LEFT";}
    if (node_col > agent_col) {move = "RIGHT";}
    return move;
};


class Explore
{
    public:
        int getRandomInt(int min, int max) 
        {
            static random_device rd;
            static mt19937 gen(rd());
            uniform_int_distribution<> dis(min, max);
            return dis(gen);
        }

        vector<string> getValidMoves(Maze* maze, Agent* agent)
        {
            vector<string> valid_moves;
            pair<int, int> a_pos = agent->getPosition();
            string maze_value = maze->getValue(a_pos.first-1, a_pos.second);
            if (maze->getValue(a_pos.first-1, a_pos.second) != "#") {valid_moves.push_back("UP");}
            if (maze->getValue(a_pos.first+1, a_pos.second) != "#") {valid_moves.push_back("DOWN");}
            if (maze->getValue(a_pos.first, a_pos.second-1) != "#") {valid_moves.push_back("LEFT");}
            if (maze->getValue(a_pos.first, a_pos.second+1) != "#") {valid_moves.push_back("RIGHT");}
            return valid_moves;
        }

        pair<int, int> updatePosition(Maze* maze, pair<int, int>* pos, string move)
        {
            pair<int, int> pos_maj = *pos;
            if (move == "UP") {pos_maj.first -= 1;}
            else if (move == "RIGHT") {pos_maj.second += 1;}
            else if (move == "DOWN") {pos_maj.first += 1;}
            else if (move == "LEFT") {pos_maj.second -= 1;}
            else 
            {
                string fname = __FUNCTION__;
                throw invalid_argument(fname + ": Invalid move " + move + " !");
            }
            return pos_maj;
        }

        bool isValidPosition(Maze* maze, pair<int, int>* pos)
        {
            bool valid_position = false;
            if (pos->first < 0 || pos->first >= maze->n_rows || pos->second < 0 || pos->second >= maze->n_cols) {return valid_position;}
            if (maze->getValue(pos->first, pos->second) != ".") {return valid_position;}
            valid_position = true;
            return valid_position;
        }

        int computeNCellsDiscovered(Maze* maze, pair<int, int>* pos)
        {
            int n_cells_discovered = 0;

            pair<int, int> s_row_up = make_pair(max(0, pos->first-2), max(0, pos->second-2));
            pair<int, int> e_row_up = make_pair(max(0, pos->first-2), min(maze->n_cols, pos->second+2));

            pair<int, int> s_row_down = make_pair(min(maze->n_rows, pos->first+2), max(0, pos->second-2));
            pair<int, int> e_row_down = make_pair(min(maze->n_rows, pos->first+2), min(maze->n_cols, pos->second+2));

            pair<int, int> s_col_left = make_pair(max(0, pos->first-2), max(0, pos->second-2));
            pair<int, int> e_col_left = make_pair(min(maze->n_rows, pos->first+2), max(0, pos->second-2));

            pair<int, int> s_col_right = make_pair(max(0, pos->first-2), min(maze->n_cols, pos->second+2));
            pair<int, int> e_col_right = make_pair(min(maze->n_rows, pos->first+2), min(maze->n_cols, pos->second+2));

            for (int i=s_row_up.second; i<=e_row_up.second; i++)
            {
                if (maze->getValue(s_row_up.first, i) == "?") {n_cells_discovered++;}
            }
            for (int i=s_row_down.second; i<=e_row_down.second; i++)
            {
                if (maze->getValue(s_row_down.first, i) == "?") {n_cells_discovered++;}
            }
            for (int i=s_col_left.first; i<=e_col_left.first; i++)
            {
                if (maze->getValue(i, s_col_left.second) == "?") {n_cells_discovered++;}
            }
            for (int i=s_col_right.first; i<=e_col_right.first; i++)
            {
                if (maze->getValue(i, s_col_right.second) == "?") {n_cells_discovered++;}
            }
            return n_cells_discovered;
        }

        string operator()(Maze* maze, Agent* agent, Timer* t)
        {
            string selected_move = "";
            int max_cells_discovered = 0;
            string moves[4] = {"UP", "DOWN", "LEFT", "RIGHT"}; 
            pair<int, int> pos = agent->getPosition();
            vector<string> valid_moves;

            // 1- search if the neighbours reveal some cells
            for (string move: moves)
            {
                pair<int, int> new_pos = updatePosition(maze, &pos, move);
                bool is_valid_position = isValidPosition(maze, &new_pos);
                if (is_valid_position)
                {
                    valid_moves.push_back(move);
                    int n_cells_discovered = computeNCellsDiscovered(maze, &new_pos);
                    if (n_cells_discovered > max_cells_discovered)
                    {
                        selected_move = move;
                        max_cells_discovered = n_cells_discovered;
                    }
                }
            }

            // 2- visit neighbours is not sufficient
            if (max_cells_discovered == 0)
            {   
                // compute difference between walkable pos and agent pos visited
                vector<pair<int, int>> relevant_positions;
                set_difference(
                    maze->walkable.begin(), 
                    maze->walkable.end(),
                    agent->visited.begin(), 
                    agent->visited.end(), 
                    back_inserter(relevant_positions)
                );

                // shuffle vector
                random_device rd;
                mt19937 g(rd());
                shuffle(relevant_positions.begin(), relevant_positions.end(), g);

                // iterate on relevant cells until out of time
                int it = 0;
                double t_exec = t->measure();
                const int LIMIT_T = 140;  // 150ms max
                while (t_exec < LIMIT_T)
                {
                    pair<int, int> r_pos = relevant_positions[it];
                    pair<int, int> target = (*maze)[r_pos.first][r_pos.second].to_pair();
                    int n_cells_discovered = computeNCellsDiscovered(maze, &target);

                    if (n_cells_discovered > 0)
                    {
                        vector<Node*> nodesPath;
                        bool reachable = bfs(
                            &(*maze)[r_pos.first][r_pos.second],
                            &(*maze)[agent->row][agent->col], 
                            maze, 
                            nodesPath
                        );
                        if (reachable)
                        {
                            selected_move = findMove(agent, nodesPath[1]);
                            max_cells_discovered = n_cells_discovered;
                            break;
                        }
                    }
                    it++;
                    t_exec = t->measure();
                }

                // 3- if no smart move has been found, take random valid move
                if (max_cells_discovered == 0)
                {
                    int rand_number = getRandomInt(0, valid_moves.size()-1);
                    selected_move = valid_moves[rand_number];
                }
            }
            return selected_move;
        }
};


int main()
{
    // env variables
    int n_rows, n_cols, a, agent_row, agent_col; 
    cin >> n_rows >> n_cols >> a; cin.ignore();

    int iteration = 0;
    int it_restart_exploring = 0;

    Timer t;
    Maze maze(n_rows, n_cols);
    Agent agent;
    Explore explore;
    string move = "";

    // game loop
    while (1) 
    {
        t.start();
        // read
        cin >> agent_row >> agent_col; cin.ignore();

        // init agent position at the first iteration only
        agent.initFirstPosition(iteration, agent_row, agent_col);

        // check if agent is on command cell
        if (maze.pos_command.first==agent.row && maze.pos_command.second==agent.col) {maze.alarm_enable = true;}

        maze.updateData(iteration);

        // check if command is reachable
        if (maze.command_found && !maze.command_reachable) 
        {
            // BFS to go to command
            pair<int, int> agent_pos = agent.getPosition();
            maze.command_reachable = isReachable(&maze[maze.pos_command.first][maze.pos_command.second], &maze[agent_pos.first][agent_pos.second], &maze);
        }

        // explore until command is reachable 
        if (!maze.command_reachable) {move = explore(&maze, &agent, &t);}

        // go to command
        else if (maze.command_reachable && !maze.alarm_enable)
        {
            // check special case: command is reachable but path between command & agent.start is too long
            vector<Node*> nodesPath;
            bfs(
                &maze[maze.pos_command.first][maze.pos_command.second], 
                &maze[agent.row_start][agent.col_start], 
                &maze, 
                nodesPath
            );
            if (nodesPath.size()>a) {move = explore(&maze, &agent, &t);}
            else
            {
                pair<int, int> agent_pos = agent.getPosition();
                bfs(
                    &maze[maze.pos_command.first][maze.pos_command.second], 
                    &maze[agent_pos.first][agent_pos.second], 
                    &maze, 
                    nodesPath
                );
                move = findMove(&agent, nodesPath[0]);
            }
        }

        // alarm is enable, come back to start point
        else if (maze.alarm_enable) 
        {
            vector<Node*> nodesPath;
            pair<int, int> agent_pos = agent.getPosition();
            bfs(
                &maze[agent.row_start][agent.col_start], 
                &maze[agent_pos.first][agent_pos.second], 
                &maze, 
                nodesPath
            );
            move = findMove(&agent, nodesPath[0]);
        }

        // real move
        agent.move(move);
        cout << move << endl; // agent's next move (UP DOWN LEFT or RIGHT).
        iteration++;
    }
}