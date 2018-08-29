for i = 1 : 24
    openfig( [ int2str( i ) '.fig' ] );
    disp( [ int2str( i ) ' ' get( gcf, 'Name' ) ] );
    set( gcf, 'WindowState', 'maximized' );
    drawnow;
    tightfig( gcf );
    set( gcf, 'WindowState', 'maximized' );
    drawnow;
    saveas( gcf, [ int2str( i ) '.emf' ] );
    close all;
end
