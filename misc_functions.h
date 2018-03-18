#ifndef misc_functions_h
#define misc_functions_h

vector<double> linspace( double start, double end, int num_of_steps );


void print( vector<double> vec );

void print( vector<vector<double> > mat );

void print( vector<vector<vector<double> > > mat );

vector<vector<vector<double> > > mesh_grid( double x_start, double x_end, int x_steps, double y_start, double y_end, int y_steps );

void hline(int n);

void dotline(int n);


void output_file(string file_name, vector<double> &vec );

void output_file(string file_name, vector< vector<double> > &mat );


//model independent string comparision
bool icompare_pred(unsigned char a, unsigned char b);
bool icompare(string const& a, string const& b);

void res_print(vector<double> &res);

double normalCDF(double value);

double normalPDF(double value);

void display_opening(string target_val, string method1, string method2 );
void display_opening(string target_val, string method1);

void display_ending();
void display_section_split();



#endif /* misc_functions_h */
