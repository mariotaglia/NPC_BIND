!#############################################################
!
! NPC program v3.14
! changes from v 3.13
! added ellipsoid
! Ptype = 4
! added uncapped cylinder
! Ptype = 5
!
! NPC program v3.13
! changes from v3.12
! read monomer properties from definitions
!
! NPC program v3.12
! changes from v3.11
! calculate the end-to-end distance of the NUPs and save them to endtoend.dat file
!
!
!
!
!
! NPC program v3.11
! changes from v3.10
! added janus spherical cargos using Ptype = 3
!
!
!
! NPC program v3.10
! changes from v3.9
! add rounded cylindrical cargoes
! add savetodiskflag = 12 to save protein properties, then stop
!
! NPC program v3.9
! changes from v3.8
! fix proteinqC
!
!
!
!
! NPC program v3.8
! changes from v3.7
!
! Allows ionizable cargo, controlled by weakP and pKaP
!

! NPC program v3.7
! changes from v3.6
!
! - SEED key in DEFINITIONS.txt 
! (edited: nanoporo.f, globals.f, DEFINITIONS.txt)
! - 
!
!
! NPC program v3.6
! changes from v3.5 
! 
! Allows using cylindrical nanoparticles by setting Ptype = 2
! (edited: nanoporo.f, globals.f, DEFINITIONS.txt)
! savetodisk_flag = 10 do not save any output
! savetodisk_flag = 11  save avpol+poten+protC
! (edited nanoporo.f



      module mprotein
      real*8, allocatable :: protein(:,:)
      real*8, allocatable :: proteinq(:,:)
      real*8, allocatable :: proteinC(:)
      real*8, allocatable :: proteinqC(:)
      real*8, allocatable :: proteinh(:, :)
      real*8, allocatable :: proteinhC(:)
      integer Ppos ! protein position in delta units
      integer PdimR ! protein radius in delta units      
      integer PdimH ! protein length in delta units, only for type 2
      integer Ptype
      real*8 prot_vol ! internal protein volume fraction
      real*8 prot_q   ! total protein charge
      real*8 prot_q2   
      real*8 hst
      real*8 hst2
      integer hrange
      real*8 pKaP, KaP, K0P
      integer weakP

      endmodule mprotein

      module posmk
      real*8, allocatable :: current(:,:)
      integer*2, allocatable :: firstcell(:,:,:)
      integer*2, allocatable :: nextbead(:)
      endmodule posmk

      module mparameters_chain
      integer N_chains ! number of chains, should be equal to size in the actual implementation
      integer, allocatable :: long(:) ! length of chain i, see chains.definitions.f
      real*8, allocatable :: chainsperdelta(:) ! number of equivalent chains grafted at position i
      real*8, allocatable :: zposition(:)      ! z position of chain i
      integer, allocatable :: segtype(:,:)
      real*8, allocatable :: endtoend(:)
      real*8 endtoend_av
      end module mparameters_chain

      module mparameters_monomer
      integer N_poorsol ! number of different kais
      integer N_monomer ! number of different monomer types
      real*8, allocatable :: st_matrix(:,:) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      integer, allocatable :: zpol(:)  ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      integer, allocatable :: hydroph(:) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      real*8, allocatable ::  pKa(:), Ka(:), K0(:)
      real*8, allocatable ::  pKbind(:), Kbind(:), rbind(:), Kbind0(:)
      real*8 pKbindread
      real*8, allocatable :: henergy(:) 
      integer monomerz
      real*8 monomerpKa
      real*8 nmonkap(100)
      character (len = 100)nmonkapbuffer
      integer, allocatable :: center(:)
      endmodule mparameters_monomer

      module mkinsol
      double precision, allocatable :: pp(:)
      endmodule 

      module mlookup
      integer*4, allocatable :: rp(:), rm(:), zp(:), zm(:)
      integer*4, allocatable :: indexa(:,:) ! 1 = R, 2 = Z
      endmodule mlookup

      module mchains
      integer*4, allocatable :: inc(:,:)
!      integer*4, allocatable :: fs(:)
!      integer*1, allocatable :: displ(:,:)
      endmodule mchains
       
      module mgraft
      real*8 posgraft(2) ! grafting point positions coord, GP
      endmodule mgraft

      module mkai
      integer Xulimit
      real*8, allocatable :: Xu(:,:,:)
      integer, allocatable :: nXu(:), Xulist_cell(:,:)
      real*8, allocatable :: Xulist_value(:,:)
      real*8 sumXu
      endmodule mkai

      module mKaps_s
      real*8, allocatable :: KapsV_s(:,:,:)
      integer, allocatable :: nKap_s(:), Kaplist_s_cell(:,:)
      real*8, allocatable :: Kaplist_s_value(:,:)
      real*8, allocatable :: v0r(:,:)
      real*8, allocatable :: v0z(:,:)
      endmodule mKaps_s

      module mKaps
      integer Kapd
      real*8, allocatable :: Kaps(:,:,:,:)
      real*8, allocatable :: KapsV(:,:,:)
      integer, allocatable :: nKap(:), Kaplist_cell(:,:)
      real*8, allocatable :: Kaplist_value(:,:,:)
      real*8 vkap
      integer, allocatable :: kapmask(:) 
      integer, parameter :: kapsrank = 0 ! processor that takes care of Kaps
      endmodule mKaps


      module mncells
      integer ncells
      endmodule mncells

      module mporesystem
      integer*4, allocatable :: matriz(:, :)
      endmodule mporesystem

      module mrands
      real*8 rands
      external rands
      integer seed
      endmodule mrands

      module msegme
      real*8, allocatable :: in1(:,:)  ! Posicion del segmento
      integer, allocatable :: celda(:)
      integer flag
      endmodule msegme

      module mvariables
! saveindex
      integer saveindex

! Volumen
      real*8 vsol
      real*8 vpol
      real*8 vsalt
! Volumen fraction
      real*8, allocatable :: avpol(:, :, :)
      real*8, allocatable :: receptor(:)
      real*8, allocatable :: funbound(:,:)
      real*8, allocatable :: fbound(:,:,:)
      real*8, allocatable :: qtot(:) 
      real*8, allocatable :: rhokap(:) 
      real*8, allocatable :: rhokapb(:) 
      real*8, allocatable :: psi2(:) 
      real*8, allocatable :: xsol(:) 
      real*8, allocatable :: xtotal2(:, :) 
      real*8, allocatable :: xpos(:) ! pos ion
      real*8, allocatable :: xneg(:) ! neg ioni
      real*8, allocatable :: xHplus(:) ! H+
      real*8, allocatable :: xOHmin(:) ! OH-
! Bulk
      real*8 xsolbulk ! volume fraction solvent in bulk
      real*8, allocatable :: xtotalbulk(:)
      real*8, allocatable :: xpotbulk(:)
      real*8, allocatable :: fdisbulk(:)

      real*8 xkapbulk
      real*8 xkapbulks(200)
      integer nxkapbulk
      real*8 xposbulk, xposbulk2, xnegbulk, xsalt,csalt
      real*8 expmupos, expmuneg
      real*8 expmukap
      real*8 chainperdelta
      real*8 sigma, sigmadelta,  mupol1, mupol2, mupold
      real*8 sigma1, sigma2
      integer*1 sigmaflag
! Charge
      real*8 zpos, zneg
      real*8 sigmaq
      real*8 lb
      real*8 constq
      real*8 betae
! Weak pol
      real*8, allocatable :: fdis(:, :)
      real*8, allocatable :: fdisP(:)
      real*8 Kw, pKw, pHbulk, expmuHplus, expmuOHmin
      real*8 xHplusbulk, cHplus, pOHbulk, xOHminbulk, cOHmin
      integer strongpoly
      real*8 shift
      real*8 norma
      integer iter
! chains
      real*8 lnpro, pro, q, sumprolnpro
! solver
      integer savetodisk_flag
! others 
      integer infile 
      endmodule mvariables

      module mparameters
      real*8 st

      integer maxlong
      integer poretype
      integer cadenastype
      integer dimR, dimZ         ! system total size is dimR x dimZ
      real*8 delta
      real*8 lseg
      integer cuantas
      integer, allocatable :: newcuantas(:)
      real*8 Na
      parameter (Na=6.02d23)    ! Avogadro's number
      real*8 pi
      parameter (pi=3.14159)    ! pi
      integer CdimR
      integer CdimZ
      integer CdimZmin
      integer BdimR
      integer TdimR
      integer RdimZ
      real*8 Curvature
      integer savetodisk_type
      integer readkai
      integer readkap
      integer save_every
      real*8 dseg
      integer nearbonds
      integer mcube
      integer ncha_max
      real*8 b2,d2                    ! here calculated once from lseg and dseg (keep unchanged)
      integer wantedCPUsecs
      integer calq
      real*8 qprob0
      integer MCsteps
      integer temp_long
      integer readseed
      endmodule mparameters


