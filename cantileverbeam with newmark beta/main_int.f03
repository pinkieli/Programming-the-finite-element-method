module main
implicit none 
interface 
SUBROUTINE deemat(dee,e,v)
!
! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
! (three dimensions).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
REAL(iwp),INTENT(IN)::e,v
 REAL(iwp),INTENT(OUT)::dee(:,:)
end subroutine deemat 

SUBROUTINE ecmat(ecm,fun,ndof,nodof)
!
! This subroutine forms the element consistent mass matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::fun(:)
 REAL(iwp),INTENT(OUT)::ecm(:,:)
 INTEGER,INTENT(IN)::nodof,ndof
end subroutine ecmat

     SUBROUTINE elmat(area,rho,emm)
!
! This subroutine forms the "analytical" lumped mass matrix for
! quadrilateral 4- or 8-node plane strain elements.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::area,rho  
 REAL(iwp),INTENT(OUT)::emm(:,:)
 end subroutine elmat 

 SUBROUTINE fkdiag(kdiag,g)
!
! This subroutine computes the skyline profile.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::g(:)
 INTEGER,INTENT(OUT)::kdiag(:)
   end subroutine fkdiag 

      SUBROUTINE formnf(nf)
!
! This subroutine forms the nf matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN OUT)::nf(:,:)
end subroutine formnf 

SUBROUTINE fsparv(kv,km,g,kdiag)
!
! This subroutine assembles element matrices into a symmetric skyline
! global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::g(:),kdiag(:)
 REAL(iwp),INTENT(IN)::km(:,:)
 REAL(iwp),INTENT(OUT)::kv(:) 
end subroutine fsparv 

 SUBROUTINE geom_rect(element,iel,x_coords,y_coords,coord,num,dir)
!
! This subroutine forms the coordinates and connectivity for a
! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
! or quadrilateral elements (4, 8 or 9-node) counting in the
! x- or y-dir. 
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:)
 REAL(iwp),INTENT(OUT)::coord(:,:)
 CHARACTER(LEN=15),INTENT(IN)::element
 CHARACTER(LEN=1),INTENT(IN)::dir
 INTEGER,INTENT(IN)::iel
 INTEGER,INTENT(OUT)::num(:)
     end subroutine geom_rect 

SUBROUTINE getname(argv,nlen)
!
! This subroutine reads the base name of data file.
!
 IMPLICIT NONE
 INTEGER::narg
 INTEGER,INTENT(OUT)::nlen
 INTEGER::lnblnk,iargc
 CHARACTER(*),INTENT(OUT)::argv
end subroutine getname

SUBROUTINE invert(matrix)
!
! This subroutine inverts a small square matrix onto itself.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::matrix(:,:)
end subroutine invert 

 SUBROUTINE linmul_sky(kv,disps,loads,kdiag)
!
! This subroutine forms the product of symmetric matrix stored as
! a skyline and a vector.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:),disps(0:)
 REAL(iwp),INTENT(OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
end subroutine linmul_sky 

SUBROUTINE mesh(g_coord,g_num,argv,nlen,ips)
!
! This subroutine produces a PostScript output file "*.msh" displaying
! the undeformed finite element mesh.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::g_coord(:,:)
 INTEGER,INTENT(IN)::g_num(:,:),ips,nlen
 CHARACTER(*),INTENT(IN)::argv
end subroutine mesh

SUBROUTINE mesh_size(element,nod,nels,nn,nxe,nye,nze)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE
 CHARACTER(LEN=15),INTENT(IN)::element
 INTEGER,INTENT(IN)::nod,nxe,nye
 INTEGER,INTENT(IN),OPTIONAL::nze
 INTEGER,INTENT(OUT)::nels,nn
end subroutine mesh_size

SUBROUTINE num_to_g(num,nf,g)
!
! This subroutine finds the g vector from num and nf.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::num(:),nf(:,:)  
 INTEGER,INTENT(OUT)::g(:)
end subroutine num_to_g

SUBROUTINE sample(element,s,wt)
!
! This subroutine returns the local coordinates and weighting coefficients
! of the integrating points.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(OUT)::s(:,:)
 REAL(iwp),INTENT(OUT),OPTIONAL::wt(:)
 CHARACTER(*),INTENT(IN)::element
end subroutine sample 

 SUBROUTINE shape_der(der,points,i)
!
!   This subroutine produces derivatives of shape functions withe respect
!   to local coordinates.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::i
 REAL(iwp),INTENT(IN)::points(:,:)
 REAL(iwp),INTENT(OUT)::der(:,:)
end subroutine shape_der

SUBROUTINE shape_fun(fun,points,i)
!
!   This subroutine computes the values of the shape functions.
!   to local coordinates
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(in)::i
 REAL(iwp),INTENT(IN)::points(:,:)
 REAL(iwp),INTENT(OUT)::fun(:)
end subroutine shape_fun

SUBROUTINE spabac(kv,loads,kdiag)
!
! This subroutine performs Cholesky forward and back-substitution
! on a symmetric skyline global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:)
 REAL(iwp),INTENT(IN OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
end subroutine spabac

SUBROUTINE sparin(kv,kdiag)
!
! This subroutine performs Cholesky factorisation on a symmetric
! skyline global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::kv(:)
 INTEGER,INTENT(IN)::kdiag(:)
end subroutine sparin

FUNCTION determinant(jac)RESULT(det)
!
! This function returns the determinant of a 1x1, 2x2 or 3x3
! Jacobian matrix.
!
 IMPLICIT NONE    
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::jac(:,:)
 REAL(iwp)::det
 end function determinant 

SUBROUTINE beemat(bee,deriv)
!
! This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::deriv(:,:)
 REAL(iwp),INTENT(OUT)::bee(:,:)
 end subroutine beemat  

end interface

end module main 

