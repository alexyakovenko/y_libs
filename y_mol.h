//This is a header of universal molecules distribution in y_system
#define Y_MOL 0x1

//Some general notes about this module:
//0. Order of elements in massives are usually has some meaning - don't shuffle objects in lists unless you know well of what are you doing.
//1. Neighbors and order are assumed to be aligned so higher-valence binders appears first in the neighbors lists.

#ifndef Y_SYSTEM
#include "y_system.h"
#endif
#ifndef Y_TXT
#include "y_txt.h"
#endif
#ifndef Y_FILE
#include "y_file.h"
#endif
#ifndef Y_LIST
#include "y_list.h"
#endif
#ifndef Y_MATH
#include "y_math.h"
#endif
#ifndef Y_VECTOR
#include "y_vector.h"
#endif
#ifndef Y_MATRIX
#include "y_matrix.h"
#endif
#ifndef Y_SMATRIX
#include "y_smatrix.h"
#endif
#ifndef Y_GRAPH
#include "y_graph.h"
#endif
#ifndef Y_GEOMETRY
#include "y_geometry.h"
#endif


//(Selected) Mendileev's table definition
#define CHEM_ATOM_TYPE_HYDROGEN     1
#define CHEM_ATOM_TYPE_CARBON       6 
#define CHEM_ATOM_TYPE_NITROGEN     7 
#define CHEM_ATOM_TYPE_OXYGEN       8 
#define CHEM_ATOM_TYPE_FLUORINE     9
#define CHEM_ATOM_TYPE_SODIUM      11  //Ionic form only 
#define CHEM_ATOM_TYPE_MAGNESIUM   12  //Ionic form only  
#define CHEM_ATOM_TYPE_SILICON     14  //Rather formal knowledge
#define CHEM_ATOM_TYPE_PHOSPHOR    15
#define CHEM_ATOM_TYPE_SULFUR      16 
#define CHEM_ATOM_TYPE_CHLORINE    17
#define CHEM_ATOM_TYPE_POTASSIUM   19  //Ionic form only
#define CHEM_ATOM_TYPE_CALCIUM     20  //Ionic form only
#define CHEM_ATOM_TYPE_BROMINE     35
#define CHEM_ATOM_TYPE_IODINE      53

//Atomic and ionic radii of the atoms
#define CHEM_ATOM_TYPE_AMOUNT      55
                                                    //e           H        He
static float CHEM_ATOM_MASS[CHEM_ATOM_TYPE_AMOUNT]={ 0.00054858, 1.008,   4.003,
//Li        Be       B        C        N        O        F        Ne
  6.940,   9.012,  10.810,  12.011,  14.007,  15.999,  18.998,  20.180,
//Na        Mg       Al       Si       P        S        Cl       Ar
 22.990,  24.305,  26.982,  28.085,  30.974,  32.060,  35.451,  39.948,
//K         Ca       Sc       Ti       V        Cr       Mn       Fe       Co       Ni       Cu       Zn       Ga       Ge       As       Se       Br       Kr
 39.098,  40.078,  44.956,  47.867,  50.942,  51.996,  54.938,  55.845,  58.933,  58.693,  63.546,  65.380,  69.723,  72.630,  74.922,  78.971,  79.904,  83.798,
//Rb        Sr       Y        Zr       Nb       Mo       Tc*      Ru       Rh       Pd       Ag       Cd       In       Sn       Sb       Tl       I        Xe
 85.468,  87.621,  88.906,  91.224,  92.906,  95.951,  98.000, 101.072, 102.906, 106.421, 107.868, 112.414, 114.818, 118.711, 121.760, 127.603, 126.904, 131.293
};
                                                     //e     H     He
static float CHEM_COV_RADIUS[CHEM_ATOM_TYPE_AMOUNT]={ 0.00, 0.31, 0.28,
//Li   Be    B     C     N     O     F     Ne
1.28, 0.96, 0.84, 0.77, 0.71, 0.66, 0.64, 0.58,
//Na   Mg    Al    Si    P     S     Cl    Ar
1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
//K    Ca    Sc   Ti     V     Cr    Mn    Fe    Co    Ni    Cu    Zn    Ga    Ge    As    Se    Br    Kr
2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.22, 1.19, 1.20, 1.20, 1.16,
//Rb   Sr    Y     Zr    Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd    In    Sn    Sb    Tl    I     Xe
2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40
};
                                                     //e     H     He
static float CHEM_VDW_RADIUS[CHEM_ATOM_TYPE_AMOUNT]={ 0.00, 1.20, 1.40,
//Li   Be    B     C     N     O     F     Ne
1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.35, 1.54,
//Na@  Mg@   Al@   Si@   P     S     Cl    Ar
2.32, 1.96, 1.75, 1.68, 1.80, 1.80, 1.75, 1.88, 
//K@   Ca@   Sc@   Ti@   V@    Cr@   Mn@   Fe@   Co@   Ni@   Cu@   Zn@   Ga@   Ge@   As@   Se    Br    Kr
2.88, 2.41, 1.98, 1.80, 1.72, 1.67, 1.66, 1.65, 1.64, 1.63, 1.73, 1.77, 1.75, 1.77, 1.81, 1.90, 1.85, 2.02,
//Rb@  Sr@   Y@    Zr@   Nb@   Mo@   Tc@   Ru@   Rh@   Pd@   Ag@   Cd@   In@   Sn@   Sb@   Tl@    I     Xe
3.04, 2.63, 2.20, 1.96, 1.86, 1.76, 1.73, 1.81, 1.75, 1.86, 1.77, 1.98, 1.96, 1.99, 2.03, 1.98, 1.98, 2.16
}; // @ - Batsanov S.S. "Van der Waals radii of elements" Inorganic Materials 37, 9, 871-885, 2001

//Bond types definition
#define CHEM_BOND_TYPE_SINGLE      '1'
#define CHEM_BOND_TYPE_DOUBLE      '2'
#define CHEM_BOND_TYPE_TRIPLE      '3'
#define CHEM_BOND_TYPE_QUATERNARY  '4'  //Formal knowledge only
#define CHEM_BOND_TYPE_AROMATIC    'r'
#define CHEM_BOND_TYPE_AMIDE       'm'


#define YSOLV_WAT_NAME "WAT"

#define MAX_ANCHOR_CYCLE_SIZE 0xA
#define MAX_ATOM_NN 0x3
#define MAX_CHEM_BOND_LENGHT 2.4

#define MAX_HB_LENGHT  3.6
#define MAX_HB_LENGHT2 MAX_HB_LENGHT*MAX_HB_LENGHT
#define MIN_HB_COS     0.5

//One residue topology representation in y_system
typedef struct{
              t_list              atoms;         //list of atoms names
              char              *ctypes;         //atoms chemical types
              unsigned int       nedges;         //number of edges in residue graph
              t_edge             *edges;         //graphs edges
              }t_res;                            //one residue topology representation structure

//Atoms force field representation in y_system
typedef struct{
              char              chem_id;         //atoms chemical type
              unsigned int        gtype;         //atoms type in gaff, NOTE the 4-th place of the name is left for envinronmental suffixes
              unsigned int        ytype;         //atoms type in yff1 
              unsigned int      generic;         //atoms generic type in yff1/gaff
              char                  ppn;         //number of pi electrons that atom donates into resonance system
              char                  pnn;         //number of nonsigma bonds
              char                   vn;         //number of the arbitrary valences
              char                   nn;         //number of coordinations
              char                   rn;         //involvement into resonance
              char                cycle;         //cycle size
              char                   nD;         //atoms inverse donors potential
              double               mass;         //atoms mass
              double               A, B;         //atoms Van-der-Waals YFF1
              double           Eng, Hrd;         //atoms Oliferenko electronegativity and hardnesses
              double               rvdw;         //Van-der-Waals atom radius
              double             sgamma;         //gamma coeficient in GB desolvation energy theory
              }t_ff_a;                           //atoms force field representation structure

//Bonds parameters
typedef struct{
              unsigned int      atom[2];
              double                k,v;
              }t_ff_b;

//Angles parameters
typedef struct{
              unsigned int      atom[3];
              double                k,v;
              }t_ff_g;

//Impropers parameters
typedef struct{
              unsigned int         type;         //type==0 means A->(B-C-D) angle, type>0 means A->(C-B-D) angle with amount around the bond to devide k-th potential
              unsigned int      atom[4];
              double               k, v;
              }t_ff_i;

//Torsions parameters
#define GAFF_GENERIC_IMPROPER_K   1.1
#define SIZE_DIH 5
typedef struct{
              unsigned int   atom[4], d;         //dihedrals atoms, number of dihedrals around the bond to devide K-th potential
              double        k[SIZE_DIH];         //dihedral potentials 
              double        v[SIZE_DIH];         //dihedral value 
              double        n[SIZE_DIH];         //dihedral harminics (up to SIZE_DIH functions per angle)
              }t_ff_d;
typedef struct{
              unsigned int      atom[4];         //dihedrals elements
              double            k[3], v;         //dihedral potentials, value and harminics 
              }t_ff_t;

//Dihedrals
typedef struct{
              unsigned int      atom[2];         //vertices of edges
              }t_ff_p;

//Topology db representation in y_system
typedef struct{
              double yff0_sA[4], yff0_sB[4], yff0_sC[4], yff0_sD[4], yff0_sQ[4],
                     yff0_delta, yff0_sF, yff0_sE; //yff0 steric and coulomb smoothing constant 
              double           **A, **B;         //yff1 vdw parameters
              unsigned int nnb_corrections;      //amount of corrected vdw parameters
              t_list              *ress;         //residues list { sizeof(unsigned int) letters per residue name is poermitted }
              t_res                *res;         //residues topologies
              unsigned int      _natoms;         //number of atoms in residues a
              char (*_atoms)[sizeof(int)];       //atoms names in residues
              char             *_ctypes;         //atoms types in residues
              unsigned int      _nedges;         //number of graph edges in residues
              t_edge            *_edges;         //graphs edges in residues
              //YFF1 parameters               
              unsigned int       size_a;         //number of yff1 atom types in topologies
              t_ff_a              *ff_a;         //yff1 atoms force field parameters
              //GAFF parameters
              unsigned int       size_n;         //amount of GAFF atom types
              unsigned int        *ff_n;         //GAFF atom names  
              unsigned int       size_b;         //number of bonds
              t_ff_b              *ff_b;         //bonds topologies
              double br[0xA][0xA], bK[0xA][0xA]; //empirical GAFF parameters for bonds
              unsigned int       size_g;         //number of angles
              t_ff_g              *ff_g;         //angles topologies
              double   gZ[0xA], gC[0xA];         //empirical GAFF parameters for angles 
              unsigned int       nitype;         //number of 1->2 impropers type
              unsigned int       size_i;         //number of known dihedral configurations
              t_ff_i              *ff_i;         //dihedral configurations
              unsigned int       size_d;         //number of known dihedral configurations
              t_ff_d              *ff_d;         //dihedral configurations GAFF
              t_ff_t              *ff_t;         //dihedral configurations
              }t_top;                            //topology db representation structure

//Molecular structure representation in y_system (compact)
//NB! A str(ucture) can't consist of no atoms, this heuristics is frequently assumed
typedef struct{
              char                    *name;     //molecule name
              t_list                  *ress;     //molecule complex sequence { sizeof(unsigne int) letters per residue name is permitted }
//            NOTE. str->ress=0x0 means that there is only one residues starting at 0 and ending at str->size_a while str->rsize keeps the residues name rather than a pointer
              unsigned int           *rsize;     //starting atom number in current residue (ress->size+1 massive lenght)
//            NOTE. rsize holds residue name for mono-residue molecules (i.e. if str->ress=0x0)
              int                 start_rid;     //the first residue in list 
              unsigned int           natoms;     //amount of atoms in the molecule          
              char   (*anames)[sizeof(int)];     //atoms names in molecule
              char                       *a;     //atoms chemical types
              t_vec                      *r;     //atoms coordinates
              unsigned int           nedges;     //number of bonds in molecule
              t_edge                 *edges;     //bonds in molecule
              }t_str;                            //edges are: if size_a>0xFE |uint|uint|char| (9 bytes) else |uchar|uchar|char| (3 bytes)
//Molecular system representation in y_system
//NB! A mol(ecule) can't consist of no atoms, this heuristics is fequently assumed
typedef struct{
//                YSTR part
              char                    *name;     //molecule name
              t_list                  *ress;     //molecule complex sequence { sizeof(unsigne int) letters per residue name is permitted }
              unsigned int           *rsize;     //starting atom number in current residue (ress->size+1 massive lenght)
              int                 start_rid;     //the first residue in list 
              unsigned int           natoms;     //amount of atoms in the molecule          
              char   (*anames)[sizeof(int)];     //atoms names in molecule
              char                       *a;     //atoms chemical types
              t_vec                      *r;     //atoms coordinates
              unsigned int           nedges;     //number of bonds in molecule
              t_edge                 *edges;     //bonds in molecule
//                Topology part
              unsigned int          size_ar;     //number of aromatic cycles in molecule
              t_clist               *cycles;     //cycles hyperlist
              t_clist              *anchors;     //lists of molecule anchors
              unsigned int          naedges;     //number of anchors edges in molecule
              t_edge                *aedges;     //anchors bond matrix
              unsigned int          *ytypes;     //atoms YFF1 types
//                Molecular Mechanics FF part
//                Electron distribution functions 
              unsigned int nvatoms, nvedges;     //number of virtual atoms and virtual bonds in the molecule
              int                   *vatoms;     //virtual atoms charges
              unsigned int     (*vedges)[2];     //virtual edges [atom_id][vatom_id]
              double           *engs, *hrds;     //atom electronegativites (real + virtual) and bond hardnesses (real + virtual) that solves the first electric moment
              char                   cmtype;     // sparse or dense triangle matrix is used for charge matrix / operator
              union     {t_smatrix *sL;
                         t_dmatrix *dL; } C;     //triangle 'charge matrix' - internal molecular operator
              double               *charges;     //atoms charges
//                FF Parameters
              unsigned int         *fftypes;     //atoms FF types [ytypes are default] (A and B are provided above the mol, at the system level) 
              unsigned int           size_b;     //Amount of bonds
              t_ff_b                  *ff_b;     //Bonds
              unsigned int           size_g;     //Amount of angles
              t_ff_g                  *ff_g;     //Angles
              unsigned int           size_i;     //Amount of impropers
              t_ff_i                  *ff_i;     //Impropers
              unsigned int           size_d;     //Amount of torsions
              t_ff_d                  *ff_d;     //Torsions
              t_list                  *excl;     //Lists of nonbonded exclusions (due to bonded interactions) per atom (partners are sorted and less than current atom_id)
              unsigned int           size_c;     //Amount of constraints
              t_ff_b                  *ff_c;     //Constraints
              }t_mol;                            //molecular system representation structure


typedef struct{
              t_edge  edge;               //Rotational edge (type is a register to travel rtree)
              unsigned int nrbranch;      //Number of anchor neighbors
              unsigned int *rbranch;      //Neighboring anchors (zero is rotational one)
              unsigned int    cross;      //Amount (from the upper side) cycles edges bonded to current anchor
              }t_rbranch;

typedef struct{
              unsigned int root;          //Root of current branch
              unsigned int  nrbranch;      //total number of branches in tree (mols anchors +1 for super-root)
              unsigned int  _nrbranch;     //number of reserved elements in _rbranch
              unsigned int *_rbranch;      //massive of in_rbranches int massives of neighbors
              t_rbranch     *rbranch;      //massive of rotational branches
              unsigned int    nidofs;      //number of active internal freeedom degreases
              }t_rtree; //NB! Only root atoms are aligned (as nonroot is not possible in general case)


//Structure record a bit modified standard str: deleted *name, *start_id, *r, (*anames) but added *q0; 
typedef struct{
              int                          id;    //molecule name
              t_list                    *ress;    //molecule complex sequence { sizeof(unsigned int) letters per residue name is permitted }
//            NOTE. dbstr->ress=0x0 means that there is only one residues starting at 0 and ending at dbstr->size_a while str->rsize keeps the residues name rather than a pointer
              unsigned int             *rsize;    //starting atom number in current residue (ress->size+1 massive lenght)
//            NOTE. rsize holds residue name for mono-residue molecules (i.e. if dbstr->ress=0x0)
              unsigned int     natoms, nedges;    //amount of vertices and edges in molecular graph
              char                         *a;    //massive of atom types (Z cores) 
              t_edge                   *edges;    //massive of edges   
              }t_dbstr;  

#ifndef Y_CHARGE
#include "y_charge.h"
#endif
#ifndef Y_MBUILD
#include "y_mbuild.h"
#endif

#define AMIDE_BOND_CONSTRAINT_K 420.

/*
 *
 * Atom types has range 0...255           [char]
 * Residue bonds number has range 0...255 [char]
 * Residue atoms number has range 0...255 [char]
 * mol has range 0...NAX_INT^2            [unsigned int]
 *
 */

/*
 * The general logic of the molecular object manipulations are:
 * 1. Upload topology as atom and residue data storage
 * 2. Upload structure in meta_data_massive and define chem types and bond orders
 * 3. Create mol object coarse by defining atom types and anchors graph
 * 4. Parameterize molecule by signing dihedral potential and charge walues
 */

#define YMOL_WAT_NAME "WAT " 

//********************************     T O P O L O G Y     P A R T     ******************************************/

//This function free top structure
void  free_top(t_top *top);

//This functions reorder members of the topological units so they gives the lowest possible value in compare_xxxx functions
//This function order atom types in the topology primitives 
inline void order_bond(unsigned int *x0,unsigned int *x1,unsigned int a0,unsigned int a1);
inline void order_angle(unsigned int *x0,unsigned int *x1,unsigned int *x2,unsigned int a0,unsigned int a1,unsigned int a2);
inline void order_impr(unsigned int type,unsigned int *x0,unsigned int *x1,unsigned int *x2,unsigned int *x3,unsigned int a0,unsigned int a1,unsigned int a2,unsigned int a3);
inline void order_tors(unsigned int *x0,unsigned int *x1,unsigned int *x2,unsigned int *x3,unsigned int a0,unsigned int a1,unsigned int a2,unsigned int a3);
//This function compare two angles topolies to decode whos is lower(greater)
int compare_bonds(const void *b1,const void *b2);
//This function compare two angles topolies to decode whos is lower(greater)
int compare_angles(const void *a1,const void *a2);
//This function compare two impropers topolies to decode who is lower(greater)
int compare_imprs(const void *i1,const void *i2);
//This function compare two torsions topolies to decode who is lower(greater)
int compare_dihs(const void *d1,const void *d2);
//This function compare two torsions topolies to decode who is lower(greater)
int compare_torss(const void *t1,const void *t2);

//This function upload topoly from files
t_top *read_top(FILE *in_atoms_yff1,FILE *in_pairs_yff1,FILE *in_atoms_gaff,FILE *in_bonds_gaff,FILE *in_angles_gaff,FILE *in_imprs_gaff,FILE *in_dihs_gaff,FILE *in_ress);

//This function performs default topology import from ../top/*
t_top *import_top();


//************************************    I / O     P R O C E S S O R S       *********************************/


//This function free memory of str
//Note it doesn't delete t_str itself
inline void free_str(register t_str *str);

//This function free mol structure
void free_mol(t_mol *mol);

//This function allocate memory for str
//Note str=0x0 means create new str
inline t_str *alloc_str(t_str *_str,register unsigned int size_r,register unsigned int size_a,register unsigned int nedges);
//This function allocate solid str
inline t_str *alloc_solid_str(register unsigned int len,register unsigned int size_r,register unsigned int size_a,register unsigned int nedges);

//This function define chemical type from the name
char name_to_chemid(register char *lexem);
//This function define chemical name from the type
int chemid_to_name(register char id);

//This function assign unique namme for atoms in str. The names should be formatted as lexems.
char get_unique_atom_name(unsigned int names_size,char (*anames)[sizeof(int)]);

//This function load molecule from mol2 file
t_str *read_mol2(FILE *in);

//This function exports mol2 files with charges
char write_mol2(FILE *out,char (*order)[4],t_mol *mol);

//This function reads first str from sdf file
//Count is the line # returned inn the case of error
t_str *read_sdf(unsigned int *count,FILE *in);

//This function import coordinates from PDBQT file (id-th structure; it return ylib_errno==YERROR_LEGAL if there is no id-th structure in the input file).
//It doesn't produce a bond matrix, just residues, atom names and corresponding coordinates
t_str *read_pdbqt(unsigned int id,FILE *in);

//This function write pdbqt file
//Note. This function requires superanchor and superrtree to save correctly 
char write_pdbqt(FILE *out,char (*order)[4],t_clist *neighbors,t_clist *sanchors,t_rtree *srtree,t_mol *mol);

//This function writes gromacs itp from parameterized mol
char export_itp(FILE *out,char (*order)[4],t_clist *neighbrs,unsigned int *anchor_id,t_mol *mol,t_top *top);

//This function writes gro file for topology generated with a function above
char write_gro(FILE *out,t_mol *mol);

//This function reads str in compact solid form. Such str should be deleted with only free(str) command;
t_str *read_solid_str(FILE *in);

//This function reads str
//NOTE. input str==0 means create new str
t_str *read_str(t_str *str,FILE *in);

//This function writes str
char write_str(FILE *out,t_str *str);

//This function reads str from a solid memory block
t_str *read_str_from_memory(void *vp);
//This function writes str into solid memory block
void *write_str_to_memory(size_t *size,t_str *str);


//This function reads ymol
t_mol *read_ymol(FILE *in);

//This function writtes ymol
char write_ymol(FILE *out,t_mol *mol);

//START OF DBSTR PART

//Tips: dbstr is solid - free(dbstr) is OK
//This function allocates solid dbstr structure
t_dbstr *alloc_solid_dbstr(unsigned int size_r,unsigned int size_a,unsigned int nedges);
//This function reads dbstr structure from file [edited from y_mol.c]
t_dbstr *read_dbstr(FILE *in);
//This function write dbstr to hdd (including the header) [edited from y_mol.c]
char write_dbstr(t_dbstr *dbstr,FILE *out);
//This function does str_to_dbstr conversion
t_dbstr *str_to_dbstr(t_str *str);
//This function does dbstr_to_str conversion (since dbstr is smaller not all data is recovered)
t_str *dbstr_to_str(t_vec  *r,t_dbstr *dbstr); 

//END OF DBSTR PART

/************************************      S T R      E D I T I O N     P A R T        *****************************************************/
//We DON'T use neighbors here becuse the molecular graph is rebuilded upon the callings 

//The order[_i] consists of 4 chars -> | amount_of_3_bonds | amount_of_2_bonds | amount of amide/quasy conjugated bonds | amount_of_1_bonds | and OPTIONALLY the neighbours in list are sorted correspondingly

//This funcion calculates chemical 'order' of atoms
char *calc_order(char (**order)[4],unsigned int size_a,unsigned int size_b,t_edge *b);

//This funcion reorder create ordered neighbors list
void arrange_neighbors(char (*order)[4],t_clist *neighbors,unsigned int nvertices,unsigned int nedges,t_edge *edges);
//This function orders neighbours list accordingly to bonds' order 
char *define_neighbors_order(t_clist *neighbors,unsigned int size_a,unsigned int nedges,t_edge *edges);

//This function is a service module for the protonate_str
//It returns the id of the atom to be inserted
inline unsigned int str_add_hydrogen(register unsigned int atom_id,register t_str *str);
//This function deletes a hydrohen from str.
//It returns the atom_id that was deleted or (unsigned int)-1 on failure 
inline unsigned int str_del_hydrogen(register unsigned int atom_id,register t_str *str);

//This is a primitive function to full-fill missing hydogens
unsigned int protonate_str(register t_str *str);

//This function edits R-(N=O)-OH -> R-N(=O)2
unsigned int update_NO2(char (*order)[4],t_str *str);

//This function edits amide bond: R1-(C-X-H)=N-R2 -> R1-(C=X)-(N-H)-R2 where X={O,S}
unsigned int update_amides(char (*order)[4],t_str *str);

//This function edits R-N=C-(NH2)2 -> R-NH-(C=NH)-NH2
//Note. It should be called after amide bond definition rutine
unsigned int update_guanidines(char (**order)[4],t_clist *_neighbors,t_str *str);

//This function edits R(0-3)-N -> R(0-3)-NH(+)
unsigned int update_amines(char (**order)[4],t_clist *_neighbors,t_str *str);

//This function edits R-CO-OH -> R-CO-O(-)
unsigned int update_carboxyles(char (*order)[4],t_str *str);

//This function edits R-SO2-OH -> R-SO2-O(-) and R-SO-OH -> R-SO-O(-)
unsigned int update_sulfaites(char (*order)[4],t_str *str);

//This function edits R-(P=O)-(OH)2 -> R->(P=O)-O2(2-)
unsigned int update_phosphates(char (*order)[4],t_str *str);

//This function edits R-(Si=O)-OH -> R->(Si=O)-O(-)
unsigned int update_silicates(char (*order)[4],t_str *str);

//
// This function need to be completed. It should resolve double bonds topology.
//
//This function correct topology of pi-electron subgraph
char update_double_bonds(char (*order)[4],t_clist *neighbors,t_str *str);


//This function do compilation of initial str (apply all the above heuristics)
char compile_str(double pH,char verbose,char (**order)[4],t_clist **neighbors,t_str *str);


//This function creates a massive of (*char)[sizeof(int)] to denote | atom_type | 0xFF ^ {n amount} | ...
unsigned int *calc_brutto_f(unsigned int natoms,char *a);

//This routine calculates brutto formula of a compound (fragment) in form of a text string
unsigned int print_brutto_f(char **brt_f,unsigned int *brutto_f);



/************************************      M O L      C O M P I L A T I O N     P A R T        *****************************************************/




//This function determines aromacity of cycles.
unsigned int classify_cycles(char (*order)[4],t_clist *neighbors,t_clist *cycles,t_str *str);

//*************   T H E     F I R S T      O R D E R       C O M P I L E R   ************/

//This function define YFF1 atom types of mol from chemical description
t_mol *compile_mol_YFF1(char (*order)[4],t_clist *neighbors,unsigned int size_ar,t_clist *cycles,t_str *str,t_top *top);

//This function implements topological proton transfering on the base of empirical rules.
//Please note that actual molecular graph should be modified previously at structure level (see numeous functions above).
//Rule 1. R-[C==[NH2]2(+)] -e is bound to each nitrogen.
//Rule 2. R3-N(sp3)-H(+) -e is bound to nitrogen.
//Rule 3. R-[O=C-O(-)] +e is bound to each oxigen.
//Rule 4. R-[O=S(+6)-O(-)] +e is bound to each oxigen.
//Rule 5. R-[O=S(+4)-O(-)] +e is bound to each oxigen
//Rule 6. R-[O=P(+5)-O(-)n] +ne are bound to each oxigen.
//Rule 7. R-[O=Si(+4)-O(-)n] +ne are bound to each oxigen.
char ionize_mol_empirically(char (*order)[4],t_clist *neighbors,t_mol *mol,t_top *top);


//*************   T H E     S E C O N D      O R D E R       C O M P I L E R S   ************/

//This function shows created anchors
void show_anchors(register t_mol *mol);
//This function defines topological part of mol - anchors and their graph. Neighbors hyper-list is used.
char disassemble_mol(unsigned int **_anchor_id,t_clist *neighbors,t_mol *mol);

//This function construct molecular mechanic force field of molecule. Neighbors and anchors hyper-lists are used.
char compose_mol(char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_adjacency *adjacency,t_clist **inacycles,t_mol *mol,t_top *top);

//GAFF parameterizators

//This function does bonds compilation
//NOTE. ff_b atoms should contain ytypes of the bond
inline char parameterize_bonds_gaff(char generic,unsigned int type,t_ff_b *ff_b,char gtypes0,char gtypes1,t_mol *mol,t_top *top);
//This function does angles compilation
#define INTERANCHOR_BONDED_MONOANGLE_HANDICAP 1.25
inline char parameterize_angles_gaff(char generic,t_ff_g *ff_g,char gtypes0,char gtypes1,char gtypes2,t_mol *mol,t_top *top);
//This function does impropers compilation
//NOTE. ff_i atoms should contain ytypes of the angle
inline char parameterize_impropers_gaff(char generic,t_ff_i *ff_i,char gtypes0,char gtypes1,char gtypes2,char gtypes3,t_mol *mol,t_top *top);
//This function does torsions compilation
//NOTE. ff_d atoms should contain ytypes of the angle
inline char parameterize_dihedrals_gaff(char generic,t_ff_d *ff_d,char gtypes0,char gtypes1,char gtypes2,char gtypes3,t_mol *mol,t_top *top);

//This function compile gaff part of the FF model
//Note. It demands that aromatic cycles can't contain X_=_Y or Y=X=Z atoms
unsigned int parameterize_mol_GAFF(char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_clist *inacycle,t_mol *mol,t_top *top);

//This function adjust parameters values to general logic, runtime problems and (optional) match given geometry
unsigned int adjust_mol_parameters(char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_clist *inacycles,t_mol *mol,t_vec *r,t_top *top);


//This function compiles a mechanical model of a molecule from the scratch (i.e. str)
t_mol *compile_mol(char verbose,double ph,char (*order)[4],t_clist *neighbors,unsigned int **_anchor_id,t_adjacency **_adjacency,t_vec *r,t_str *str,t_top *top);





// ******************************        M O L E C U L E S    A N D    S U B S T R U C T U R E     C O N S T R U C T I O N S     P A R T    ******************************  
//NB!! Only single bond can separate two anchors
//NB!! Each gag starts with root atom, the edges are [acceptor_id],[root_id]:[1]



//This function construct water molecule t_mol structure
t_str *construct_wat(t_str *str,t_top *top);


//Gags trigger works as follows:
// | geometry bit | topology bit | output bit* |
// |       1/0    |       1/0    |     1/0     |
// *Means 0 - output atoms number, 1 - output bonds number

//This is service function to define approximate single bond lenght between two atom with defined chemical types
//The supported types H,C,N,O,F,Si,P,S,Cl,Br,I; all other pairs are very approximately detemined on their row in periodic table
double get_approximate_single_bond_lenght(char chem_type_1,char chem_type_2);

// C A R B O N    G A G S
//This function makes carboxyl gag insted of an atom0->atom1 :: atom0->COO(-)
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-C = 1.55342A, C-O = 1.23445A, X(C)-C-O = 115.3Deg, O-C-O = 129.4Deg
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.5268)
//O1  0.0000  1.4044  1.6639
//O2  0.0000 -1.4044  1.6639 
void make_carboxyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes cyane gag insted of an atom0->atom1 :: atom0->[NH2-C(+)=NH2]
//It adds 7 atoms and 7 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.49999A, C-N = 1.3085A, X(C)-C-N = 120.0Deg, N-H = 0.99803A, C-N-H = 121.6Deg, all in one plane
//X    0.0000  0.00000  0.00000
//C    0.0000  0.00000  1.00000 (1.49999)
//N1   0.0000  1.13319  1.65425
//H11  0.0000  1.16106  2.65189
//H12  0.0000  2.01111  2.12894
//N2   0.0000 -1.13319  1.65425
//H21  0.0000 -1.16106  2.65189
//H22  0.0000 -2.01111  2.12894
void make_cguanidine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes cyane gag insted of an atom0->atom1 :: atom0->C_=_N
//It adds 2 atoms and 2 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all on OZ line, X(C)-C = 1.467A, C-N = 1.13478A, X(C)-C-N = 180.0Deg
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.5268)
//N   0.0000  0.0000  2.1348
void make_cyane_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes amide gag insted of an atom0->atom1 :: atom0->CO-NH2
//It adds 5 atoms and 5 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz-plane, X(C)-C = 1.51314A, C-O = 1.19775A, X(C)-C-O = 122.8Deg, C-N = 1.35592A, X(C)-C-N = 115.0Deg, O-C-N = 122.2Deg, N-H = 0.9920A, C-N-H = 120.4Deg, H-N-H = 119.2Deg 
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.51314)
//O   0.0000  1.0068  1.6488 
//N   0.0000 -1.2105  1.5645
//H1  0.0000 -0.3548  2.0664 
//H2  0.0000 -2.0664  2.0664
void make_camide_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes methyl gag insted of an atom0->atom1 :: atom0->CH3
//It adds 4 atoms and 4 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.52679A, C-H = 1.08579A, X(C)-C-H = 111.217Deg, z-dih = PI*2/3, H-C-H = 107.671Deg
//X   0.0000  0.0000  0.0000
//C   0.0000  0.0000  1.0000 (1.5268)
//H1  0.0000  1.0122  1.3930
//H2 -0.8766 -0.5061  1.3930 
//H3  0.8766 -0.5061  1.3930
void make_methyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes ethene gag insted of an atom0->atom1 :: R->CH=CH2
//It adds 5 atoms and 5 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-C = 1.51213A, C-H = 1.08782A, X(C)-C-H = 116.6Deg, X(C)-C-C = 124.8Deg, C-C = 1.31853A, C-H = 1.0775A, C-C-H = 121.6Deg, all in the same plane
//X    0.0000  0.00000  0.00000
//C1   0.0000  0.00000  1.00000 (1.51213)
//H1   0.0000 -0.97268  1.48708
//C2   0.0000  1.08271  1.75250
//H21  0.0000  1.02256  2.82832 
//H22  0.0000  2.07009  2.18388
void make_ethene_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);

// N I T R O G E N    G A G S
//This function makes charged amine gag insted of an atom0->atom1 :: atom0->NH3(+)
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: CH3-NH3 = 1.50585A, N-H = 1.01021A, X(CH3)-N-H = 111.472Deg, z-dih = PI*2/3, H-N-H = 107.398Deg
//X   0.0000   0.0000  0.0000
//N   0.0000   0.0000  1.0000 (1.50585)
//H1  0.0000   0.9401  1.3698
//H2  0.81415 -0.47005 1.3698 
//H3 -0.81415 -0.47005 1.3698 
void make_camine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes charged nguanidine gag insted of an atom0->atom1 :: atom0->NH~(HN=C(+)~NH2)
//It adds 9 atoms and 9 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-NH = 1.46284A, N-H = 0.996A, X(CH3)-N-H = 117.40Deg, N-C = 1.31841A, C(X)-N-C = 125.17Deg, C-N = 1.32676A, N-C-N = 120.0Deg, C-N-H = 121.60Deg, all in one plane.
//X    0.00000  0.00000  0.00000
//N0   0.00000  0.00000  1.00000 (1.46284)
//H0   0.00000 -0.88426  1.45836
//C    0.00000  1.07773  1.75941 
//N1   0.00000  0.95817  3.08077 
//H11  0.00000  0.06628  3.52410 
//H12  0.00000  1.75601  2.48456 
//N2   0.00000  2.28184  2.31655
//H21  0.00000  3.11172  2.86730 
//H22  0.00000  2.39925  3.30561 
void make_nguanidine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes amide gag insted of an atom0->atom1 :: atom0->NH2
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(CO)-NH2 = 1.34798A, N-H = 0.993921A, X(CO)-N-H = 119.053Deg, H-N-H = 119.336Deg
//X   0.0000  0.0000  0.0000
//N   0.0000  0.0000  1.0000 (1.3480)
//H1  0.0000  0.8689  1.4827
//H2  0.0000 -0.8689  1.4827 
void make_pamine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes amide gag insted of an atom0->atom1 :: atom0->N=CH2
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-NH2 = 1.45565A, N-C = 1.2809A, X(C)-N-C = 122.47Deg, C-H = 1.0875A, N-C-H = 121.65Deg, all in one plane
//X   0.0000  0.00000  0.00000
//N   0.0000  0.00000  1.00000 (1.45565)
//C   0.0000  1.08030  1.68823
//H1  0.0000  1.06417  2.77561
//H2  0.0000  2.05898  2.16240 
void make_nimine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes namide gag insted of an atom0->atom1 :: atom0->NH-(C=O)H
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-N = 1.44772A, N-H = 0.9931A, X(C)-N-H = 119.56Deg, X(C)-N-C = 121.76Deg, H-N-C = 118.68Deg, N-C = 1.34474A, C=O = 1.19569A, C-H = 1.09128A, N-C-O = 124.76Deg, N-C-H = 113.0Deg, O-C-H = 122.24Deg, H-N-C-O = 180.0Deg
//X   0.0000  0.00000  0.0000
//N   0.0000  0.00000  1.0000  (1.44772)
//H   0.0000  0.86384  1.48993
//C   0.0000 -1.14338  1.70782
//O   0.0000 -2.12570  2.38953
//H   0.0000 -1.13885  2.13422 
void make_namide_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes nitro gag insted of an atom0->atom1 :: atom0-> NO2
//It adds 3 atoms and 3 edges to str
//Note. r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, all symmetrical on OZ, X(C)-N = 1.47864A, N-O = 1.19220A, X(C)-N-O = 117.1Deg, O-N-O = 125.8Deg
//X   0.0000  0.00000  0.0000
//N   0.0000  0.00000  1.0000  (1.47864)
//O   0.0000  1.06131  1.5431
//O   0.0000 -1.06131  1.5431
void make_nitro_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);

// O X Y G E N    G A G S 
//This function makes hydroxyl gag insted of an atom0->atom1 :: atom0->OH
//It adds 2 atoms and 2 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-O = 1.39965A, O-H = 0.94630A, X(C)-O-H = 109.45Deg
//X   0.0000  0.0000  0.0000
//O   0.0000  0.0000  1.0000 (1.39965)
//H   0.0000  0.8923  1.3151
void make_hydroxyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes methoxy gag insted of an atom0->atom1 :: atom0->O-CH3
//It adds 2 atoms and 2 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all but H2&H3 in yz plane, X(C)-O = 1.3922A, O-C = 1.3922A, X(C)-O-C = 115.70Deg, C-H = 1.0849A, O-C-H = 109.7Deg, O->C__L_H=120Deg.
//X   0.00000  0.00000  0.00000
//O   0.00000  0.00000  1.00000 (1.3922)
//C   0.00000  1.25448  1.60374
//H1  0.00000  2.02696  0.84197
//H2  0.20185  2.17857  1.91914
//H3 -0.20185  2.17857  1.91914
void make_methoxy_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);

// S I L I C O N     G A G S 
//This function makes silane gag insted of an atom0->atom1 :: atom0->Si-H3
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-Si = 1.88831, Si-H = 1.08615, X(C)-Si-H = 111.12Deg, H-Si-H = 107.78Deg, H...X(C)-Si..H = 2*PI/3
//X    0.00000  0.00000  0.00000
//Si   0.00000  0.00000  1.00000 (1.88831)
//H1   0.00000  1.01319  1.39136
//H2  -0.87745 -0.50660  1.39136
//H3   0.87745 -0.50660  1.39136
void make_silane_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes silyl gag insted of an atom0->atom1 :: atom0->Si-(CH3)3
//It adds 13 atoms and 13 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-Si = 1.89364, Si-C = 1.89364, X(C)-Si-C = 109.47Deg, C-H = 1.08748A, Si-C-H = 111.46Deg, H-C-H = 107.41Deg, C...X(C)-Si..C = 2*PI/3, H...Si-C...H = 2*PI/3
//X     0.00000  0.00000  0.00000
//Si    0.00000  0.00000  1.00000 (1.89364)
//C1    0.00000  1.78535  1.63118
//C2   -1.54616 -0.89268  1.63118
//C3    1.54616 -0.89268  1.63118
//H11   0.00000  2.87232  1.59778
//H12  -0.28649  2.70693  2.13242
//H13   0.28649 -2.70693  2.13242
//H21  -2.48750 -1.43616  1.59778
//H22  -2.20102 -1.60157  2.13242
//H23  -2.48750 -1.10535  2.13242
//H31   2.48750 -1.43616  1.59778
//H32   2.48750 -1.10535  2.13242
//H33   2.20102 -1.60157  2.13242
void make_silyl_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);

// P H O S P H O R     G A G S 
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->P-H2
//It adds 3 atoms and 3 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.86035, P-H = 1.40416, X(C)-P-H = 98.65Deg, H-P-H = 95.09Deg, H...X(C)-P..H = 96.5Deg
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.86035)
//H1   0.00000  1.38819  1.21118
//H2  -1.37926 -0.15715  1.21118
void make_hphosphine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->P-(CH3)2
//It adds 9 atoms and 9 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.85344, P-C = 1.85344, X(C)-P-C = 100.03Deg, C-P-C = 100.03Deg, C...X(C)-P..C = 102.2Deg, C-H = 1.08534A, P-C-H = 109.8Deg, H-C-H = 107.67Deg, H...P->C-H = 2PI/3
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.85344)
//C1   0.00000  1.82606  1.32297
//C2  -1.78482 -0.38589  1.32297  
//H11  0.00000  2.36594  0.38143
//H12 -1.63996  3.42353  2.12402
//H13  1.63996  3.42353  2.12402
//H21 -2.31250 -0.49998  0.38143
//H22 -2.99964 -2.32640  2.12402
//H23 -3.69277  0.87945  2.12402
void make_mphosphine_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes phosphane gag insted of an atom0->atom1 :: atom0->P=CH2
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.8624A, P-C = 1.64612A, X(C)-P-C = 103.825Deg, C-H = 1.07666A, P-C-H = 122.17Deg, H-C-H = 115.66Deg, H..(P-C)..H = PI.
//X   0.00000  0.00000  0.00000
//P   0.00000  0.00000  1.00000 (1.89757)
//C   0.00000  1.59843  1.39335
//H1  0.00000  1.93730  2.41529
//H2  0.00000  2.37285  2.14133
void make_phosphane_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->(P=O)-O(-)2
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.89757A, P-O = 1.5190A, X(C)-P-O = 103.15Deg, O...X(C)-P..O = 2PI/3
//X   0.00000  0.00000  0.00000
//P   0.00000  0.00000  1.00000 (1.89757)
//O1  0.00000  1.47917  1.34557
//O2 -1.28100 -0.73958  1.34557
//O3  1.28100 -0.73958  1.34557
void make_phosphate2_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->(O=P~O(-))-CH3
//It adds 7 atoms and 7 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.84644A, X(C)-P-C = 101.32Deg, C-H = 1.087324A, P-C-H = 111.60Deg, P-O = 1.48887A, X(C)-P-O = 107.88Deg, C/O...X(C)-P..O = 2PI/3
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.81915)
//C1   0.00000  1.81052  1.36243
//H11  0.00000  2.00456  2.43230
//H12  0.87552  2.30222  0.94535 
//H13 -0.87552  2.30222  0.94535
//O1   1.22712 -0.70848  1.45707 
//O2  -1.22712 -0.70848  1.45707
void make_phosphate1_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes phosphate gag insted of an atom0->atom1 :: atom0->(P=O)-(CH3)2
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-P = 1.81915A, X(C)-P-C = 104.80Deg, C-H = 1.08524A, P-C-H = 109.44Deg, H..P-C..H = 2PI/3, P-O = 1.47371A, X(C)-P-O = 113.80Deg, C...X(C)-P..O = 2PI/3
//X    0.00000  0.00000  0.00000
//P    0.00000  0.00000  1.00000 (1.81915)
//C1   0.00000  1.75880  1.48513
//H11  0.00000  1.92393  2.56681
//H12  0.88627  2.31606  1.08268 
//H13 -0.88627  2.31606  1.08268
//C2   1.52316 -0.87940  1.48513 
//H21  1.66617  0.96197  2.56681
//H22  1.56263 -1.92556  1.08268
//H23  2.44890  0.39050  1.08268 
//O   -1.16774 -0.67419  1.59471
void make_phosphate0_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);

// S U L F U R     G A G S 
//This function makes sulfide gag insted of an atom0->atom1 :: atom0->S-CH3
//It adds 5 atoms and 5 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.8087A, S-C = 1.8087A, X(C)-S-C = 100.0Deg, C-H = 1.08304A, S-C-H = 111.16Deg, H-C-H = 109.5Deg, H...S-C..H = 2PI/3
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.8087)
//C   0.00000  1.78122  1.31408
//H1  0.00000  2.31454  0.37145
//H2 -1.59915  3.26166  2.09345
//H3  1.59915  3.26166  2.09345
void make_sulfide_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes thiol gag insted of an atom0->atom1 :: atom0->S-H
//It adds 2 atoms and 2 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: all in yz plane, X(C)-S = 1.81768A, S-H = 1.32663A, X(C)-S-H = 97.9Deg
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.81768)
//H   0.00000  1.31404  1.18234
void make_thiol_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes sulfite gag insted of an atom0->atom1 :: atom0->(S=O)-O(-)
//It adds 3 atoms and 3 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.82014A, S-O = 1.4944A, X(C)-S-O = 101.0Deg, O-S-O = 113.82Deg, O...X(C)-S..O = 117.1
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.82014)
//O1  0.00000  1.46694  1.28514
//O2 -1.30589 -0.66826  1.28514
void make_isulfite_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes methyl sulfite gag insted of an atom0->atom1 :: atom0->(S=O)-CH3
//It adds 6 atoms and 6 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.79647A, C-H = 1.08300A, S-C-H = 109.68Deg, X(C)-S-C = 97.66Deg, S-O = 1.48524, X(C)-S-O = 106.72Deg, C...X(C)-S..O = 110.9Deg
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.82014)
//C   0.00000  1.78044  1.23946
//H1  0.00000  2.00598  2.29871
//H2  0.88312  2.20987  0.78275
//H3 -0.88312  2.20987  0.78275
//O2 -1.32886 -0.50744  1.42730
void make_msulfite_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes sulfate gag insted of an atom0->atom1 :: atom0->(O=S=O)-O(-)
//It adds 4 atoms and 4 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.78722A, S-O = 1.45473A, X(C)-S-O = 104.46Deg, O...X(C)-S..O = 2PI/3
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.78722)
//O1  0.00000  1.40865  1.36325
//O2 -1.21992 -0.70432  1.36325
//O3  1.21992 -0.70432  1.36325
void make_isulfate_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);
//This function makes methyl-sulfate gag insted of an atom0->atom1 :: atom0->(O=S=O)-CH3
//It adds 7 atoms and 7 edges to str
//Note.  r is direction of the bond, acceptor_id is the atom to accept the group
//The geometry: X(C)-S = 1.78722A, X(C)-S-C = 104.32Deg, C-H = 1.08124A, S-C-H = 109.76Deg, H-(S->C)-H = 2PI/3, S-O = 1.45473A, X(C)-S-O = 104.46Deg, O...X(C)-S..O = 2PI/3
//X   0.00000  0.00000  0.00000
//S   0.00000  0.00000  1.00000 (1.78722)
//C1  0.00000  1.73060  1.44628
//H1  0.00000  1.83048  2.52289
//H2  0.88124  2.21162  1.04489
//H3 -0.88124  2.21162  1.04489
//O2 -1.21992 -0.70432  1.36325
//O3  1.21992 -0.70432  1.36325
void make_msulfate_gag(char mode,unsigned int *ress_id,unsigned int *atom_id,unsigned int *edge_id,t_vec *r,unsigned int acceptor_id,t_str *str);

//The gags are named as separate residues  
//The gags are one of the list:
//  function_name           chemistry             gag_name   description
//  make_carboxy_gag()    : R<-[O=C~O(-)]          [cxyg]  (carboxylic group is kept)
//  make_cguanidine_gag() : R<-[NH2~C(+)=NH2]      [cgug]  (guanidine group is kept)
//  make_cyane_gag()      : R<-C_=_N               [cyng]  (cyane group is kept)
//  make_camide_gag()     : R<-(C=O)~NH2           [cmdg]  (C near/in resonance)
//  make_methyl_gag()     : R<-CH3                 [metg]  (generic aliphatic C)
//  make_ethene_gag()     : R<-CH=CH2              [ethg]  (alkene carbon) 
//  make_camine_gag()     : R<-NH3+                [camg]  (aliphatic charged N)
//  make_nguanidine_gag() : R<-XN~[NH2~C(+)=NH2]   [ngug]  (guanidine group is kept)
//  make_pamine_gag()     : R<-NH2                 [pamg]  (N near aromatic if X is sp2 atom)
//  make_nimine_gag()     : R<-N=CH2               [pimg]  (generic sp2 N)
//  make_namide_gag()     : R<-NH~(C=O)H           [nadg]  (N near resonance if X is sp3 atom)
//  make_nitro_gag()      : R<-NO2                 [no2g]  (NO2 group is kept)
//  make_hydroxy_gag()    : R<-OH                  [hxyg]  (hydroxy gag) 
//  make_methoxy_gag()    : R<-O-CH3               [omeg]  (generic oxygen gag)
//  make_silane_gag()     : R<-SiH3                [slag]  (silane is kept)
//  make_silyl_gag()      : R<-Si(CH3)3            [slyg]  (silyls are generic silicone gag)
//  make_hphosphine_gag() : R<-PH2                 [hphg]  (3-coordinated phosphor bound to hydrogen)
//  make_mphosphine_gag() : R<-P-(CH3)2            [mphg]  (3-coordinated phosphor bound to methyls)
//  make_phosphane_gag()  : R<-P=CH2               [ppng]  (2-coordinated P near double bond)
//  make_phosphate2_gag() : R<-(P=O)~[O(-)]2       [pp2g]  (5-coordinated phosphor ion -2) 
//  make_phosphate1_gag() : R<-(H3C-P=O)~O(-)      [pp1g]  (5-coordinated phosphor ion -1) 
//  make_phosphate0_gag() : R<-(O=P)-(CH3)2        [pp0g]  (5-coordinated phosphor) 
//  make_sulfide_gag()    : R<-S-CH3               [sidg]  (2-coordinated sulfur)
//  make_thiol_gag()      : R<-SH                  [thig]  (2-coordinated sulfur bound to hydrogen)
//  make_isulfite_gag()   : R<-(S=O)~O(-)          [isig]  (4-coordinated sulfur ion)
//  make_msulfite_gag()   : R<-(S=O)-CH3           [msig]  (4-coordinated sulfur)
//  make_isulfate_gag()   : R<-(O=S=O)~O(-)        [isag]  (6-coordinated sulfur ion)
//  make_msulfate_gag()   : R<-(O=S=O)-CH3         [msag]  (6-coordinated sulfur)
//The gags bits modes (4+2+1 is a valid mode) 
// 1 - info
// 2 - topology
// 4 - geometry
//
//NB!! Only single (nonresonance, nonaromatic, nonamide etc) bond is allowed to separate anchors
//This function calls specific gag depending on it's name and pass the parameters through
char use_gags(unsigned int gag_name,char mode,unsigned int *res_name,unsigned int *atom_num,unsigned int *edge_num,t_vec *r,unsigned int acceptor_id,t_str *str);

//This function identifies the most appropriate gag for _i -> _j bond
//Note. It returns TRUE on success and FALSE it he generic gag was used
char _identify_gag(unsigned int *gag_name,unsigned int *atom_num,unsigned int *edge_num, register unsigned int _id, register unsigned int _i, register unsigned int _j,int *charged,char (*order)[4],t_clist *neighbors,t_mol *mol);

//This function split molecule into a collection of mono- and bi-root strs
char splitt_anchors(unsigned int *mrr_size,t_str ***mrr_str,unsigned int *brr_size,t_str ***brr_str,unsigned int maxerrcount,char (*order)[4],t_clist *neighbors,unsigned int *anchor_id,t_adjacency *adjacency,t_mol *mol,unsigned int default_gag_name);


/*************************************   R T R E E     P A R T    **********************************/





//This function alloc rotation tree data structure
t_rtree *alloc_rtree(unsigned int nrbranch,unsigned int nnrbranch);

//This function writes rtree to hdd 
char write_rtree(FILE *out,t_rtree *rtree);

//This function reads rtree from hdd
t_rtree *read_rtree(FILE *in);

//This function build rtree for given anchors complex list and massive of atomic edges
t_rtree *build_rtree(unsigned int root,unsigned int *anchor_id,t_clist *anchors,unsigned int nedges,t_edge *edges);

//This function changes the root of a rtree
char change_root(unsigned int new_root,t_rtree *rtree);

//This function builds super-anchors stucture from given rtree and anchors list
t_clist *build_sanchors(unsigned int size_a,unsigned int *anchor_id,t_clist *anchors,t_rtree *rtree);

//This function computes supertree (i.e. joins anchors so there is no crosses)
t_rtree *build_srtree(unsigned int sroot,unsigned int size_a,t_clist *sanchors,unsigned int nedges,t_edge *edges,t_rtree *rtree);

//This function create superroot for given rtree and active_a list.
//It cleans the active list from cycles and crosses.
//char compile_superroot(unsigned int *nrbranch,unsigned int **rbranch,t_list *active_a,t_rtree *rtree);



/*************************************   M I S S C E L L A N E O U S    P A R T    **********************************/


//This function gather an achors of given residues list
t_list *get_ress_anchors(t_list *ress,t_mol *mol);


