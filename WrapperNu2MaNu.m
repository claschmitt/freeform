function MatLabNurbs = WrapperNu2MaNu(nurbs)

if isequal(nurbs.form,'B-Spline-Surface')
    MatLabNurbs.form         = 128;
    MatLabNurbs.dim(1,1)     = nurbs.dimU;
    MatLabNurbs.dim(1,2)     = nurbs.dimV;
    MatLabNurbs.number(1,1)  = nurbs.numberU;
    MatLabNurbs.number(1,2)  = nurbs.numberV;
    MatLabNurbs.coefs(1,:,:) = reshape(nurbs.coefs(1,:),nurbs.numberU,nurbs.numberV);
    MatLabNurbs.coefs(2,:,:) = reshape(nurbs.coefs(2,:),nurbs.numberU,nurbs.numberV);
    MatLabNurbs.coefs(3,:,:) = reshape(nurbs.coefs(3,:),nurbs.numberU,nurbs.numberV);
    MatLabNurbs.coefs(4,:,:) = ones(nurbs.numberU,nurbs.numberV);
    MatLabNurbs.order(1,1)   = nurbs.orderU;
    MatLabNurbs.order(1,2)   = nurbs.orderV;
    MatLabNurbs.knots{1}     = nurbs.knotsU;
    MatLabNurbs.knots{2}     = nurbs.knotsV;
else
end