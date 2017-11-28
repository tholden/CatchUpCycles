function IsZero = AssertIsZero( In )
    IsZero = abs( In ) < 1e-12;
    if ~IsZero
        disp( In );
    end
    assert( IsZero );
end
