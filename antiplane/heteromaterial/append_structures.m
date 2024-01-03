function joinedstruct = append_structures(struct1,struct2)

joinedstruct = [];
joinedstruct.N = struct1.N+struct2.N;
joinedstruct.x2c = [struct1.x2c;struct2.x2c];
joinedstruct.x3c = [struct1.x3c;struct2.x3c];
joinedstruct.x2 = [struct1.x2;struct2.x2];
joinedstruct.x3 = [struct1.x3;struct2.x3];
joinedstruct.W = [struct1.W;struct2.W];
joinedstruct.dip = [struct1.dip;struct2.dip];
joinedstruct.nvec = [struct1.nvec;struct2.nvec];

end