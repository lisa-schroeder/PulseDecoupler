
"generates basis operators dependent on the number of spins\n
output: iu,ix,iy,iz,ip,im,ia,ib"
function basis(nspins)

  # %product operator set of basis matrices for nspins spin 1/2 nuclei
  # %nspins has to be specified in main program

  # %Complex definition
  i=sqrt(Complex(-1));

  # %Matrices for single spin 1/2
  mix=0.5*[0 1; 1 0];
  miy=0.5*[0 -i; i 0];
  miz=0.5*[1 0; 0 -1];
  mip=[0 1;0 0];
  mim=[0 0;1 0];
  mia=[1 0;0 0];
  mib=[0 0;0 1];
  ione=[1 0;0 1];

  # %Initializations
  norder=2^nspins;
  global iu=zeros(ComplexF64,norder,norder,nspins);
  global ix=zeros(ComplexF64,norder,norder,nspins);
  global iy=zeros(ComplexF64,norder,norder,nspins);
  global iz=zeros(ComplexF64,norder,norder,nspins);
  global ip=zeros(ComplexF64,norder,norder,nspins);
  global im=zeros(ComplexF64,norder,norder,nspins);
  global ia=zeros(ComplexF64,norder,norder,nspins);
  global ib=zeros(ComplexF64,norder,norder,nspins);

  # %Calculation individual spin basis matrices
  for ispins=1:nspins
    links=diagm(ones(2^(ispins-1)))
    rechts=diagm(ones(2^(nspins-ispins)))
    iu[:,:,ispins]=kron(links,kron(ione,rechts))
    ix[:,:,ispins]=kron(links,kron(mix,rechts))
    iy[:,:,ispins]=kron(links,kron(miy,rechts))
    iz[:,:,ispins]=kron(links,kron(miz,rechts))
    ip[:,:,ispins]=kron(links,kron(mip,rechts))
    im[:,:,ispins]=kron(links,kron(mim,rechts))
    ia[:,:,ispins]=kron(links,kron(mia,rechts))
    ib[:,:,ispins]=kron(links,kron(mib,rechts))
  end

  return iu,ix,iy,iz,ip,im,ia,ib
end
