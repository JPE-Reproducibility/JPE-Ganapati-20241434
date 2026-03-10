function [ ] = Figure_Sequencer( d, o ,title,suffix )


AA = figure;
p_overlay_elasticity( d,o,[  'rho'       suffix]  ,['../../Output/' title '_'] );
p_overlay_elasticity( d,o,[  'epsilon'   suffix]  ,['../../Output/' title '_'] );
p_overlay_elasticity( d,o,[  'extensive' suffix]  ,['../../Output/' title '_'] );
p_overlay_elasticity( d,o,[  'intensive' suffix]  ,['../../Output/' title '_'] );
p_overlay_elasticity( d,o,[  'theta'     suffix]  ,['../../Output/' title '_'] );
close(AA);

