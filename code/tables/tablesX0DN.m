// line 27 requires downloading the appropriate cmf data

intrinsic X0DNdata(DBound::RngIntElt, NBound::RngIntElt) -> Any
    {Outputs data table for the shimura curves X0(D;N), where D a quaternion discriminant and N coprime to D in the box D =< DBound and N <= NBound.}

    filename := Sprintf("./data/genera-tables/SignatureTableX0DN_%o_%o.txt", DBound, NBound);
    fprintf filename, "label|Glabel|all_degree1_points_known|autmuO_norms|bad_primes|cm_discriminants|coarse_class|coarse_class_num|coarse_index|coarse_label|coarse_num|conductor|curve_label|deg_mu|dims|discB|discO|fine_label|fine_num|fuchsian_index|galEnd|generators|genus|genus_minus_rank|gerbiness|aut_gerbiness|has_obstruction|index|is_coarse|is_split|lattice_labels|lattice_x|level|level_is_squarefree|level_radical|level_is_prime|level_is_prime_power|log_conductor|models|mu_label|mults|name|newforms|nu2|nu3|nu4|nu6|num_bad_primes|num_known_degree1_noncm_points|num_known_degree1_points|obstructions|order_label|parents|parents_conj|pointless|power|psl2label|q_gonality|q_gonality_bounds|qbar_gonality|qbar_gonality_bounds|ram_data_elts|rank|reductions|scalar_label|simple|squarefree|torsion|trace_hash|traces\n";
    fprintf filename, "text|text|boolean|integer[]|integer[]|integer[]|text|integer|integer|text|integer|integer[]|text|integer|integer[]|integer|integer|text|integer|integer|text|integer[]|integer|integer|integer|integer|smallint|integer|boolean|boolean|text[]|integer[]|integer|boolean|integer|boolean|boolean|numeric|smallint|text|integer[]|text|text[]|integer|integer|integer|integer|integer|integer|integer|integer[]|text|text[]|integer[]|boolean|boolean|text|integer|integer[]|integer|integer[]|numeric[]|integer|text[]|text|boolean|boolean|integer[]|bigint|integer[]\n\n";
    // Preload AutmuO sizes from precomputed polarized data for the N=1 case (deg_mu=1 entries).
    pol_autmuO_size := AssociativeArray(Integers());
    pol_mu          := AssociativeArray(Integers());
    pol_file := Open("./data/quaternion-orders/quaternion-orders-polarized.txt", "r");
    _ := Gets(pol_file); // header row
    _ := Gets(pol_file); // type row
    _ := Gets(pol_file); // blank line
    pol_line := Gets(pol_file);
    while not IsEof(pol_line) do
        if #pol_line gt 0 then
            fields := Split(pol_line, "|");
            if #fields ge 6 and StringToInteger(fields[4]) eq 1 then
                D_key := StringToInteger(fields[2]);
                pol_autmuO_size[D_key] := StringToInteger(fields[6]);
                pol_mu[D_key]          := fields[3];
            end if;
        end if;
        pol_line := Gets(pol_file);
    end while;
    delete pol_file;

    // we only want D square-free, with an even number of prime factors.
    for D in [D : D in [6..DBound] | MoebiusMu(D) eq 1] do
        // we want N that are coprime to D and square-free.
        for N in [N : N in [1..NBound] | GCD(D,N) eq 1] do
            X := CMFLoad(D, N : cmfdatafile := "./code/jacobian_decomp/cmfdata.txt", levelbound:=D*N);
            O := QuaternionOrder(D,N);
            Gamma := FuchsianGroup(O : VerifyEichler := false);
            // e := EllipticInvariants(Gamma);
            Glabel:="1.1";
            all_degree1_points_known:="\\N";
            autmuO_norms:="\\N";
            bad_primes:="\\N";
            cm_discriminants:="\\N";
            coarse_class:="a";
            coarse_class_num:="1";
            if N eq 1 and IsDefined(pol_autmuO_size, D) then
                size_AutmuO := pol_autmuO_size[D];
            else
                _, mu := HasPolarizedElementOfDegree(O,1);
                AutmuO:=Domain(Aut(O,mu));
                size_AutmuO := #AutmuO;
            end if;
            coarse_index:=size_AutmuO;
            genus:=Genus(Gamma);
            newforms, mults, dims, conductor, rank, simple := JLDecomposition(D, N, X : g := genus);
            squarefree := &and[m eq 1 : m in mults] select "T" else "F";
            newforms := #newforms eq 0 select "{}" else "{" cat Join(newforms, ",") cat "}";
            mults := #mults eq 0 select "{}" else "{" cat Join([Sprintf("%o", m) : m in mults], ",") cat "}";
            dims := #dims eq 0 select "{}" else "{" cat Join([Sprintf("%o", d) : d in dims], ",") cat "}";
            conductor := #conductor eq 0 select "{}" else "{" cat Join(["{" cat Sprintf("%o", fac[1]) cat "," cat Sprintf("%o", fac[2]) cat "}" : fac in conductor], ",") cat "}";
            coarse_label:=Sprintf("1.%o.%o.%o.1", size_AutmuO, genus, coarse_class);
            coarse_num:=1;
            curve_label:="\\N";
            deg_mu:=1;
            discB:=D;
            discO:=D*N;
            fine_label:=coarse_label;
            fine_num:="\\N";
            fuchsian_index:=size_AutmuO; // Double check that this is the correct index
            galEnd:="\\N";
            generators:="\\N";
            genus_minus_rank:=genus - rank;
            if N eq 1 and IsDefined(pol_mu, D) then
                mu_coords := eval ReplaceString(ReplaceString(pol_mu[D], "{", "["), "}", "]");
                mu_elem := O!mu_coords;
                K_seq, _ := SemidirectToNormalizerKernel(O, mu_elem);
                gerbiness := #K_seq;
                aut_gerbiness := gerbiness/2;
            else
                gerbiness := "\\N";
                aut_gerbiness := "\\N";
            end if;
            has_obstruction:="\\N";
            index:=size_AutmuO;
            is_coarse:="T";
            is_split:="\\N";
            lattice_labels:="\\N";
            lattice_x:="\\N";
            level:=N;
            level_is_squarefree:=IsSquarefree(N) select "T" else "F";
            level_radical:=&*PrimeDivisors(N);
            level_is_prime := IsPrime(N) select "T" else "F";
            level_is_prime_power := (N gt 1 and IsPrimePower(N)) select "T" else "F";
            log_conductor:="\\N";
            models:="\\N";
            name:=N eq 1 select Sprintf("X(%o;1)",D) else "X_0(" cat Sprintf("%o;%o)", D,N);
            nu2:=EllipticInvariants(Gamma)[1][2];
            nu3:=EllipticInvariants(Gamma)[2][2];
            nu4:=0;
            nu6:=0;
            num_bad_primes:=#PrimeDivisors(D*N);
            num_known_degree1_noncm_points:="\\N";
            num_known_degree1_points:="\\N";
            obstructions:="\\N";
            if N eq 1 then
                order_label:=Sprintf("%o",D);
            else
                order_label:=Sprintf("%o.%o",D,N);
            end if;
            mu_label:=order_label cat ".1";
            label:=mu_label cat "." cat coarse_label;
            parents:="\\N";
            parents_conj:="\\N";
            pointless:="\\N";
            power:="\\N";
            psl2label:=label;
            gonality_temp:=GonalityBoundListX0DN(D,N);
            q_gonality:=gonality_temp[1];
            q_gonality_bounds:=Sprintf("{%o,%o}",gonality_temp[2][1],gonality_temp[2][2]);
            qbar_gonality:=gonality_temp[3];
            qbar_gonality_bounds:=Sprintf("{%o,%o}",gonality_temp[4][1],gonality_temp[4][2]);
            ram_data_elts:="\\N";
            reductions:="\\N";
            scalar_label:="\\N";
            torsion:="\\N";
            trace_hash:="\\N";
            traces:="\\N";
            fprintf filename, Sprintf("%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o\n",label,Glabel,all_degree1_points_known,autmuO_norms,bad_primes,cm_discriminants,coarse_class,coarse_class_num,coarse_index,coarse_label,coarse_num,conductor,curve_label,deg_mu,dims,discB,discO,fine_label,fine_num,fuchsian_index,galEnd,generators,genus,genus_minus_rank,gerbiness,aut_gerbiness,has_obstruction,index,is_coarse,is_split,lattice_labels,lattice_x,level,level_is_squarefree,level_radical,level_is_prime,level_is_prime_power,log_conductor,models,mu_label,mults,name,newforms,nu2,nu3,nu4,nu6,num_bad_primes,num_known_degree1_noncm_points,num_known_degree1_points,obstructions,order_label,parents,parents_conj,pointless,power,psl2label,q_gonality,q_gonality_bounds,qbar_gonality,qbar_gonality_bounds,ram_data_elts,rank,reductions,scalar_label,simple,squarefree,torsion,trace_hash,traces);
        end for;
    end for;
    return Sprint("Table produced :)");
end intrinsic;