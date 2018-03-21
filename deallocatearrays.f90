subroutine deallocarrays

use params
deallocate(x)
deallocate(v)
deallocate(Vmat)
deallocate(Amat)
deallocate(Bmat)
deallocate(BmatInv)
deallocate(Tmat)
deallocate(Hmat)
deallocate(Psi)
deallocate(Prob)
deallocate(IPIV)
deallocate(WORK1)
deallocate(W)
deallocate(Z)
deallocate(AB)
deallocate(WORK2)
deallocate(IWORK)

end subroutine deallocarrays
