$# LS-DYNA Keyword file created by LS-PrePost 4.2 (Beta) - 23Mar2015(09:00)
$# Created on May-16-2017 (17:39:17)
*KEYWORD
*TITLE
$#                                                                         title
LS-DYNA keyword deck by LS-PrePost
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Control Cards		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
*CONTROL_CPU
$   CPUTIM
$#  cputim      
     0.000
*CONTROL_ENERGY
$     HEGN      RWEN    SLNTEN     RYLEN
$#    hgen      rwen    slnten     rylen     
         2         2         2         1
*CONTROL_HOURGLASS
$      IHQ        QH
$#     ihq        qh  
         10.10000000
*CONTROL_OUTPUT
$    NPOPT    NEECHO    NREFUP    IACCOP     OPIFS    IPNINT    IKEDIT    IFLUSH
$#   npopt    neecho    nrefup    iaccop     opifs    ipnint    ikedit    iflush
         1         3         0         0     0.000         0         0         0
$    IPRTF
$#   iprtf    ierode     tet10    msgmax    ipcurv      gmdt   ip1dblt      eocs
         0         0         2        50         0     0.000         0         0
$#   tolev    newleg    frfreq     minfo    msgflg      
         2         0         1         0         0
*CONTROL_PARALLEL
$     NCPU    NUMRHS     CONST      PARA
$#    ncpu    numrhs     const      para    
         1         0         2         0
*CONTROL_SHELL
$   WRPANG     ESORT     IRNXX    ISTUPD    THEORY       BWC     MITER      PROJ
$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj
     0.000         1         0         0         2         2         1         0
$  ROTASCL    INTGRD    LAMSHT
$# rotascl    intgrd    lamsht    cstyp6    tshell      
     0.000         0         0         1         0
$# psstupd   sidt4tu     cntco    itsflg    irquad      
         0         0         0         0         2
$#  nfail1    nfail4   psnfail    keepcs     delfr   drcpsid    drcprm      
         0         0         0         0         0         0 1.0000000
*CONTROL_TERMINATION
$   ENDTIM    ENDCYC     DTMIN    ENDENG    ENDMAS
$#  endtim    endcyc     dtmin    endeng    endmas      
0.20000000         0     0.000     0.000     0.000
*CONTROL_TIMESTEP
$   DTINIT    TSSFAC      ISDO    TSLIMT     DT2MS      LCTM     ERODE     MSIST
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
     0.0000.89999998         0     0.000     0.000         0         0         0
$   DT2MSF
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     
     0.000         0         0                         0.000
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Databases		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
*DATABASE_ABSTAT
$      dt
$#      dt    binary      lcur     ioopt     
0.00100000         0         0         1
*DATABASE_GLSTAT
$      dt
$#      dt    binary      lcur     ioopt     
0.00100000         0         0         1
*DATABASE_MATSUM
$      dt
$#      dt    binary      lcur     ioopt     
0.00100000         0         0         1
*DATABASE_RWFORC
$      dt
$#      dt    binary      lcur     ioopt     
0.00100000         0         0         1
*DATABASE_BINARY_D3PLOT
$       dt      lcdt      beam     npltc    psetid    istats    tstart      iavg
$#      dt      lcdt      beam     npltc    psetid      
0.00500000         0         0         0         0
$#   ioopt     
         0
*DATABASE_BINARY_D3THDT
$       dt      lcdt      beam     npltc    psetid    istats    tstart      iavg
$#      dt      lcdt      beam     npltc    psetid      
 1.0000000         0         0         0         0
*DATABASE_BINARY_INTFOR
$       dt      lcdt      beam     npltc    psetid    istats    tstart      iavg
$#      dt      lcdt      beam     npltc    psetid      
 1.0000000         0         0         0         0
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf       vid     death     birth
         4         3         2         1 1.0000000         01.00000E28     0.000
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         1         0         0         0         1         0         0         0
*SET_NODE_LIST_TITLE
NODESET(SPC) 1
$#     sid       da1       da2       da3       da4    solver      
         1     0.000     0.000     0.000     0.000MECH
$#    nid1      nid2      nid3      nid4      nid5      nid6      nid7      nid8
        13         2        14         3        21         1        16         4
        15         0         0         0         0         0         0         0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         2         0         0         1         0         0         0         0
*SET_NODE_LIST_TITLE
NODESET(SPC) 2
$#     sid       da1       da2       da3       da4    solver      
         2     0.000     0.000     0.000     0.000MECH
$#    nid1      nid2      nid3      nid4      nid5      nid6      nid7      nid8
        19        29        18        10        23        11        14         2
        13         0         0         0         0         0         0         0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         3         0         1         0         0         0         0         0
*SET_NODE_LIST_TITLE
NODESET(SPC) 3
$#     sid       da1       da2       da3       da4    solver      
         3     0.000     0.000     0.000     0.000MECH
$#    nid1      nid2      nid3      nid4      nid5      nid6      nid7      nid8
        19        11        14         3        24         6        20        12
        15         0         0         0         0         0         0         0
*SET_NODE_LIST
$#     sid       da1       da2       da3       da4    solver      
         4     0.000     0.000     0.000     0.000MECH
$#    nid1      nid2      nid3      nid4      nid5      nid6      nid7      nid8
         8        17        33        26        18         7        19         6
        20         0         0         0         0         0         0         0
*SET_PART_LIST
$#     sid       da1       da2       da3       da4    solver      
         1     0.000     0.000     0.000     0.000MECH
$#    pid1      pid2      pid3      pid4      pid5      pid6      pid7      pid8
         1         0         0         0         0         0         0         0
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Parts		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
*PART
$name
$#                                                                         title
Matrix
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1         1         0         0         0         0         0
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Sections		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
*SECTION_SOLID
$       id   elformu       aet
$#   secid    elform       aet   
         1         1         4
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Materials		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
*MAT_ELASTIC
$      MID        RO         E        PR        DA        DB         K
$#     mid        ro         e        pr        da        db  not used        
         1     0.000     0.000     0.000     0.000     0.000         0
*PART
$#                                                                         title
Fiber
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2         2         1         0         0         0         0         0
*SECTION_BEAM
$#   secid    elform      shrf   qr/irid       cst     scoor       nsm   
         2         3 1.0000000         2         0     0.000     0.000
$#       a     rampt    stress      
7.85398E-7     0.000     0.000
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp      
         1         0 1.0000000 1.0000000     0.000     0.000         0
$#                a1                  o1  
               0.000               0.000
          0.10000000          0.25000000
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Elements		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
*ELEMENT_SOLID (ten nodes format)
$#   eid     pid   
      57       1
$#    n1      n2      n3      n4      n5      n6      n7      n8      n9     n10
      17       8      26       5       9      22      27      25       0       0
      58       1
      18       7      26       8      10      23      27      22       0       0
      59       1
      15       4      25      12       3      21      27      24       0       0
      60       1
      20       5      26       6      12      25      27      24       0       0
      61       1
      16       9      25       4       1      22      27      21       0       0
      62       1
      13       2      23      10       1      21      27      22       0       0
      63       1
      19       6      26       7      11      24      27      23       0       0
      64       1
      14       3      24      11       2      21      27      23       0       0
*ELEMENT_BEAM
$#   eid     pid      n1      n2      n3     rt1     rr1     rt2     rr2   local
      65       2      28      29       0       0       0       0       0       0
      66       2      30      31       0       0       0       0       0       0
      67       2      32      33       0       0       0       0       0       0
*NODE
$   Node               x               y               z      tc      rc
$#   nid               x               y               z      tc      rc  
       1      1.00000000   -5.960460e-08     -1.00000000       0       0
       2   -5.960460e-08     -1.00000000     -1.00000000       0       0
       3     -0.99999899    2.980230e-07     -1.00000000       0       0
       4    1.788140e-07      1.00000000     -1.00000000       0       0
       5    2.384190e-07      1.00000000      1.00000000       0       0
       6     -1.00000000    1.788140e-07      1.00000000       0       0
       7   -4.768370e-07     -1.00000000      1.00000000       0       0
       8      0.99999899   -5.364420e-07      1.00000000       0       0
       9      1.00000000      1.00000000           0.000       0       0
      10      0.99999899     -1.00000000           0.000       0       0
      11     -1.00000000     -1.00000000           0.000       0       0
      12     -0.99999899      1.00000000           0.000       0       0
      13      1.00000000     -1.00000000     -1.00000000       0       0
      14     -1.00000000     -1.00000000     -1.00000000       0       0
      15     -0.99999899      1.00000000     -1.00000000       0       0
      16      1.00000000      1.00000000     -1.00000000       0       0
      17      1.00000000      0.99999899      1.00000000       0       0
      18      0.99999899     -1.00000000      1.00000000       0       0
      19     -1.00000000     -1.00000000      1.00000000       0       0
      20     -0.99999899      1.00000000      1.00000000       0       0
      21    5.960460e-08    1.192090e-07     -1.00000000       0       0
      22      0.99999899   -2.682210e-07    6.039610e-14       0       0
      23   -2.831220e-07     -1.00000000   -1.465500e-14       0       0
      24     -1.00000000    2.533200e-07   -6.661320e-15       0       0
      25    2.235170e-07      1.00000000   -1.065810e-14       0       0
      26   -1.192090e-07   -1.639130e-07      1.00000000       0       0
      27   -3.631357e-07   -9.934200e-09    9.492407e-15       0       0
      28   -5.960460e-08     -1.00000000     -1.00000000       0       0
      29   -4.768370e-07     -1.00000000      1.00000000       0       0
      30    5.960460e-08    1.192090e-07     -1.00000000       0       0
      31   -1.192090e-07   -1.639130e-07      1.00000000       0       0
      32    1.788140e-07      1.00000000     -1.00000000       0       0
      33    2.384190e-07      1.00000000      1.00000000       0       0
$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$		Writing Boundary Conditions		
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$
$BOUNDARY CONDITIONS
*END
