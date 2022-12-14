## map atomistic structure to Martini v2.0 beads for umodified Cellulose fibers##

package require mol
mol load pdb 7CNF.pdb

# Set bead types B1 and B3 after finding COM

foreach seg {C1	C2	C3	C4	C5	C6	C7	C8	C9	C10	C11	C12	C13	C14	C15	C16	C17	C18	C19	C20	C21	C22	C23	C24	C25	C26	C27	C28	C29	C30	C31	C32	C33	C34	C35	C36	C37	C38	C39	C40	C41	C42	C43	C44	C45	C46	C47	C48	C49	C50	C51	C52	C53	C54	C55	C56	C57	C58	C59	C60	C61	C62	C63	C64	C65	C66	C67	C68	C69	C70	C71	C72	C73	C74	C75	C76	C77	C78	C79	C80	C81	C82	C83	C84	C85	C86	C87	C88	C89	C90	C91	C92	C93	C94	C95	C96	C97	C98	C99	C100	C101	C102	C103	C104	C105	C106	C107	C108	C109	C110	C111	C112	C113	C114	C115	C116	C117	C118	C119	C120	C121	C122	C123	C124	C125	C126} {
 for {set res 1} {$res < 21} {incr res 1} {
 set b3 [atomselect top "(name O6 C6 C5) and segname $seg and resid $res"]
 set com1 [measure center $b3 weight mass]
 set c6 [atomselect top "name C6 and segname $seg and resid $res"]
 $c6 moveto $com1
 $c6 set name B3
 set b1 [atomselect top "(name C2 C3 O2 O3) and segname $seg and resid $res"]
 set com2 [measure center $b1 weight mass]
 set o2 [atomselect top "name O2 and segname $seg and resid $res"]
 $o2 moveto $com2
 $o2 set name B1
 }
}
# Set bead type B2 for all residues except 1 and 20 in each chain after finding COM

foreach seg {C1	C2	C3	C4	C5	C6	C7	C8	C9	C10	C11	C12	C13	C14	C15	C16	C17	C18	C19	C20	C21	C22	C23	C24	C25	C26	C27	C28	C29	C30	C31	C32	C33	C34	C35	C36	C37	C38	C39	C40	C41	C42	C43	C44	C45	C46	C47	C48	C49	C50	C51	C52	C53	C54	C55	C56	C57	C58	C59	C60	C61	C62	C63	C64	C65	C66	C67	C68	C69	C70	C71	C72	C73	C74	C75	C76	C77	C78	C79	C80	C81	C82	C83	C84	C85	C86	C87	C88	C89	C90	C91	C92	C93	C94	C95	C96	C97	C98	C99	C100	C101	C102	C103	C104	C105	C106	C107	C108	C109	C110	C111	C112	C113	C114	C115	C116	C117	C118	C119	C120	C121	C122	C123	C124	C125	C126} {
 for {set res 2} {$res < 20} {incr res 1} {
 set b2 [atomselect top "(name C4 O5 C1 O4) and segname $seg and resid $res"]
 set com1 [measure center $b2 weight mass]
 set c4 [atomselect top "name C4 and segname $seg and resid $res"]
 $c4 moveto $com1
 $c4 set name B2
 }
}

# Set bead type B2 for residue 1 in each chain after finding COM

foreach seg {C1	C2	C3	C4	C5	C6	C7	C8	C9	C10	C11	C12	C13	C14	C15	C16	C17	C18	C19	C20	C21	C22	C23	C24	C25	C26	C27	C28	C29	C30	C31	C32	C33	C34	C35	C36	C37	C38	C39	C40	C41	C42	C43	C44	C45	C46	C47	C48	C49	C50	C51	C52	C53	C54	C55	C56	C57	C58	C59	C60	C61	C62	C63	C64	C65	C66	C67	C68	C69	C70	C71	C72	C73	C74	C75	C76	C77	C78	C79	C80	C81	C82	C83	C84	C85	C86	C87	C88	C89	C90	C91	C92	C93	C94	C95	C96	C97	C98	C99	C100	C101	C102	C103	C104	C105	C106	C107	C108	C109	C110	C111	C112	C113	C114	C115	C116	C117	C118	C119	C120	C121	C122	C123	C124	C125	C126} {
 for {set res 1} {$res < 2} {incr res 1} {
 set b2 [atomselect top "(name O4 C4 O5 C1 O4) and segname $seg and resid $res"]
 set com1 [measure center $b2 weight mass]
 set c4 [atomselect top "name C4 and segname $seg and resid $res"]
 $c4 moveto $com1
 $c4 set name B2
 }
}

# Set bead type B2 for residue 20 in each chain after finding COM

foreach seg {C1	C2	C3	C4	C5	C6	C7	C8	C9	C10	C11	C12	C13	C14	C15	C16	C17	C18	C19	C20	C21	C22	C23	C24	C25	C26	C27	C28	C29	C30	C31	C32	C33	C34	C35	C36	C37	C38	C39	C40	C41	C42	C43	C44	C45	C46	C47	C48	C49	C50	C51	C52	C53	C54	C55	C56	C57	C58	C59	C60	C61	C62	C63	C64	C65	C66	C67	C68	C69	C70	C71	C72	C73	C74	C75	C76	C77	C78	C79	C80	C81	C82	C83	C84	C85	C86	C87	C88	C89	C90	C91	C92	C93	C94	C95	C96	C97	C98	C99	C100	C101	C102	C103	C104	C105	C106	C107	C108	C109	C110	C111	C112	C113	C114	C115	C116	C117	C118	C119	C120	C121	C122	C123	C124	C125	C126} {
 for {set res 20} {$res < 21} {incr res 1} {
 set b2 [atomselect top "(name C4 O5 C1 O1) and segname $seg and resid $res"]
 set com1 [measure center $b2 weight mass]
 set c4 [atomselect top "name C4 and segname $seg and resid $res"]
 $c4 moveto $com1
 $c4 set name B2
 }
}

set towrite [atomselect top "name B1 B2 B3"]
$towrite writepdb test.pdb
$towrite writepsf test.psf

mol load psf test.psf pdb test.pdb

set existing_B3type [atomselect top "name B3"]
$existing_B3type set name P1
	
set existing_B1type [atomselect top "name B1"]
$existing_B1type set name P4

set existing_B2type [atomselect top "name B2"]
$existing_B2type set name PX

set sel_all [atomselect top all]
$sel_all set resname CEL1X

$sel_all writepdb base-CNF-bundle.pdb
$sel_all writegro base-CNF-bundle.gro

exit
