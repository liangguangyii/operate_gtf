module type_defination
    implicit none
    type atom_type
        character(len=2):: atom_name
        integer:: atom_index
        real*8:: x, y, z, charge
    end type atom_type
    
    type gtf_type
        integer:: gtf, atom
        real*8:: exponent
    end type gtf_type
end module type_defination
    
module variable_defination
    implicit none
    integer:: s2f(-5:5, 21)=0
    real*8:: expcutoff=-40D0
    real*8:: conv6d5d(6,5), conv10f7f(10,7), conv15g9g(15,9), conv21h11h(21,11)
    integer:: shell_generacy(-5:5) = (/ 11,9,7,5,4,1,3,6,10,15,21 /) 
    integer:: type2ix(56)=(/ 0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5 /)
    integer:: type2iy(56)=(/ 0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0, 0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
    integer:: type2iz(56)=(/ 0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0, 5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
    character(len=2)::periodic_table(0:120)=(/ "Bq","H ","He", &   !Bq(number 0) is ghost atom. Bq is recorded in .fch, but X is not recorded
                    "Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
                    "Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
                    "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
                    "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
                    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
                    "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
                    "Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
                    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og","Un","Ux" /)
    real*8,parameter :: pi=3.141592653589793D0
    contains
    
    !NOTE: initialize the the variables that are different in different formats
    !init_fch for fch format
    subroutine init_fch()
        !s2f: the index of generate shell for higher augular momentum
        s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
        s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
        s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
        s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
        s2f(-1,1:4)=(/ 1,2,3,4 /)
        s2f(0,1)=1
        s2f(1,1:3)=(/ 2,3,4 /)
        s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
        s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /)
        s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
        s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)
        
        !conv****: the transform matrix, from spherical repr. to Car. repr.
        conv6d5d=0D0
        conv10f7f=0D0
        conv15g9g=0D0
        conv21h11h=0D0
        
        !conv6d5d
        ! From 5D: D 0,D+1,D-1,D+2,D-2
        ! To 6D:  1  2  3  4  5  6
        !        XX,YY,ZZ,XY,XZ,YZ
        conv6d5d(1:3,1)=(/ -0.5D0, -0.5D0, 1D0 /)
        conv6d5d(5,2)=1D0
        conv6d5d(6,3)=1D0
        conv6d5d(1:2,4)=(/ sqrt(3D0)/2D0, -sqrt(3D0)/2D0 /)
        conv6d5d(4,5)=1D0
        
        !conv10f7f
        ! From 7F: F 0,F+1,F-1,F+2,F-2,F+3,F-3
        ! To 10F:  1   2   3   4   5   6   7   8   9  10      
        !         XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ (Gaussian sequence, not identical to Multiwfn)
        conv10f7f(3,1)=1D0
        conv10f7f(6,1)=-1.5D0/sqrt(5D0)
        conv10f7f(9,1)=-1.5D0/sqrt(5D0)
        conv10f7f(1,2)=-sqrt(3D0/8D0)
        conv10f7f(4,2)=-sqrt(3D0/40D0)
        conv10f7f(7,2)=sqrt(6D0/5D0)
        conv10f7f(2,3)=-sqrt(3D0/8D0)
        conv10f7f(5,3)=-sqrt(3D0/40D0)
        conv10f7f(8,3)=sqrt(6D0/5D0)
        conv10f7f(6,4)=sqrt(3D0)/2D0
        conv10f7f(9,4)=-sqrt(3D0)/2D0
        conv10f7f(10,5)=1D0
        conv10f7f(1,6)=sqrt(5D0/8D0)
        conv10f7f(4,6)=-3D0/sqrt(8D0)
        conv10f7f(2,7)=-sqrt(5D0/8D0)
        conv10f7f(5,7)=3D0/sqrt(8D0)
        
        !conv15g9g
        ! From 9G: G 0,G+1,G-1,G+2,G-2,G+3,G-3,G+4,G-4
        ! To 15G:   1    2    3    4    5    6    7    8
        !         ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ
        !           9   10   11   12   13   14   15
        !         XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
        conv15g9g(1,1)=1D0
	    conv15g9g(3,1)=-3D0*sqrt(3D0/35D0)
	    conv15g9g(5,1)=3D0/8D0
	    conv15g9g(10,1)=-3D0*sqrt(3D0/35D0)
	    conv15g9g(12,1)=3D0/4D0*sqrt(3D0/35D0)
	    conv15g9g(15,1)=3D0/8D0
	    conv15g9g(6,2)=2D0*sqrt(5D0/14D0)
	    conv15g9g(8,2)=-1.5D0/sqrt(14D0)
	    conv15g9g(13,2)=-1.5D0*sqrt(5D0/14D0)
	    conv15g9g(2,3)=2D0*sqrt(5D0/14D0)
	    conv15g9g(4,3)=-1.5D0*sqrt(5D0/14D0)
	    conv15g9g(11,3)=-1.5D0/sqrt(14D0)
	    conv15g9g(3,4)=-3D0*sqrt(3D0/28D0)
	    conv15g9g(5,4)=sqrt(5D0)/4D0
	    conv15g9g(10,4)=3D0*sqrt(3D0/28D0)
	    conv15g9g(15,4)=-sqrt(5D0)/4D0
	    conv15g9g(7,5)=3D0/sqrt(7D0)
	    conv15g9g(9,5)=-sqrt(5D0/28D0)
	    conv15g9g(14,5)=-sqrt(5D0/28D0)
	    conv15g9g(8,6)=-3D0/sqrt(8D0)
	    conv15g9g(13,6)=sqrt(5D0/8D0)
	    conv15g9g(4,7)=-sqrt(5D0/8D0)
	    conv15g9g(11,7)=3D0/sqrt(8D0)
	    conv15g9g(5,8)=sqrt(35D0)/8D0
	    conv15g9g(12,8)=-3D0/4D0*sqrt(3D0)
	    conv15g9g(15,8)=sqrt(35D0)/8D0
	    conv15g9g(9,9)=-sqrt(5D0)/2D0
	    conv15g9g(14,9)=sqrt(5D0)/2D0
        
        !conv21h11h
        ! From 11H: H 0,H+1,H-1,H+2,H-2,H+3,H-3,H+4,H-4,H+5,H-5
        ! To 21H:   1     2     3     4     5     6     7     8     9    10
        !         ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ 
        !          11    12    13    14    15    16    17    18    19    20    21
        !         XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
        conv21h11h(1,1)=1D0
        conv21h11h(12,1)=-5D0/sqrt(21D0)
        conv21h11h(3,1)=-5D0/sqrt(21D0)
        conv21h11h(19,1)=5D0/8D0
        conv21h11h(5,1)=5D0/8D0
        conv21h11h(14,1)=sqrt(15D0/7D0)/4D0
        conv21h11h(7,2)=sqrt(5D0/3D0)
        conv21h11h(16,2)=-3D0*sqrt(5D0/28D0)
        conv21h11h(9,2)=-3D0/sqrt(28D0)
        conv21h11h(21,2)=sqrt(15D0)/8D0
        conv21h11h(11,2)=sqrt(5D0/3D0)/8D0
        conv21h11h(18,2)=sqrt(5D0/7D0)/4D0
        conv21h11h(2,3)=sqrt(5D0/3D0)
        conv21h11h(4,3)=-3D0*sqrt(5D0/28D0)
        conv21h11h(13,3)=-3D0/sqrt(28D0)
        conv21h11h(6,3)=sqrt(15D0)/8D0
        conv21h11h(20,3)=sqrt(5D0/3D0)/8D0
        conv21h11h(15,3)=sqrt(5D0/7D0)/4D0
        conv21h11h(12,4)=sqrt(5D0)/2D0
        conv21h11h(3,4)=-sqrt(5D0)/2D0
        conv21h11h(19,4)=-sqrt(35D0/3D0)/4D0
        conv21h11h(5,4)=sqrt(35D0/3D0)/4D0
        conv21h11h(8,5)=sqrt(5D0/3D0)
        conv21h11h(17,5)=-sqrt(5D0/12D0)
        conv21h11h(10,5)=-sqrt(5D0/12D0)
        conv21h11h(16,6)=sqrt(5D0/6D0)
        conv21h11h(9,6)=-sqrt(1.5D0)
        conv21h11h(21,6)=-sqrt(17.5D0)/8D0
        conv21h11h(11,6)=sqrt(17.5D0)/8D0
        conv21h11h(18,6)=sqrt(5D0/6D0)/4D0
        conv21h11h(4,7)=-sqrt(5D0/6D0)
        conv21h11h(13,7)=sqrt(1.5D0)
        conv21h11h(20,7)=-sqrt(17.5D0)/8D0
        conv21h11h(6,7)=sqrt(17.5D0)/8D0
        conv21h11h(15,7)=-sqrt(5D0/6D0)/4D0
        conv21h11h(19,8)=sqrt(35D0)/8D0
        conv21h11h(5,8)=sqrt(35D0)/8D0
        conv21h11h(14,8)=-0.75D0*sqrt(3D0)
        conv21h11h(17,9)=sqrt(5D0)/2D0
        conv21h11h(10,9)=-sqrt(5D0)/2D0
        conv21h11h(21,10)=3D0/8D0*sqrt(3.5D0)
        conv21h11h(11,10)=5D0/8D0*sqrt(3.5D0)
        conv21h11h(18,10)=-1.25D0*sqrt(1.5D0)
        conv21h11h(6,11)=3D0/8D0*sqrt(3.5D0)
        conv21h11h(20,11)=5D0/8D0*sqrt(3.5D0)
        conv21h11h(15,11)=-1.25D0*sqrt(1.5D0)
    end subroutine init_fch
end module variable_defination    