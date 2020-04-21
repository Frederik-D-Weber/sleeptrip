function st = sleepStage2hypnNum(st,unknownIsNaN,plot_yaxequidist)
switch st
    case {'W' 'w' 'Wake' 'wake' 'WAKE' '0'}
        st = 0;
    case {'S1' 'N1' 'stage 1' 'Stage 1' 'Stage1' 'STAGE 1' 'STAGE1' 'S 1' 'Stadium 1' 'STADIUM 1' 'STADIUM1' '1'}
        if plot_yaxequidist
            st = -2;
        else
            st = -1;
        end
    case {'S2' 'N2' 'stage 2' 'Stage 2' 'Stage2' 'STAGE 2' 'STAGE2' 'S 2' 'Stadium 2' 'STADIUM 2' 'STADIUM2' '2' }
        if plot_yaxequidist
            st = -3;
        else
            st = -2;
        end
    case {'S3' 'N3' 'stage 3' 'Stage 3' 'Stage3' 'STAGE 3' 'STAGE3' 'S 3' 'Stadium 3' 'STADIUM 3' 'STADIUM3' '3' 'SWS'}
        if plot_yaxequidist
            st = -4;
        else
            st = -3;
        end
    case {'S4' 'N4' 'stage 4' 'Stage 4' 'Stage4' 'STAGE 4' 'STAGE4' 'S 4' 'Stadium 4' 'STADIUM 4' 'STADIUM4' '4' 'SWS4'}
        if plot_yaxequidist
            st = -5;
        else
            st = -4;
        end
    case {'REM' 'R' 'r' 'Rem' 'rem' 'Stage 5' 'Stage5' 'STAGE 5' 'STAGE5' 'S 5' 'Stadium 5' 'STADIUM 5' 'STADIUM5' '5'}
        if plot_yaxequidist
            st = -1;
        else
            st = -0.5;
        end
    case {'MT' 'mt' 'movement' 'Movement' 'Movement Time' 'MovementTime' '8'}
        if plot_yaxequidist
            st = 1;
        else
            st = 0.5;
        end
    case {'A' 'a' 'Artifact' 'Artefact' 'artifact' 'artefact' 'Artf' 'Artif.'}
        if plot_yaxequidist
            st = 2;
        else
            st = 1;
        end
    case {'?' '???' 'unknown', 'Unknown' '-'  'X' 'x'}
        if unknownIsNaN
            st = NaN;
        else
            if plot_yaxequidist
                st = 3;
            else
                st = 1.5;
            end
        end
    otherwise
        if unknownIsNaN
            st = NaN;
        else
            if plot_yaxequidist
                st = 3;
            else
                st = 1.5;
            end;
        end
end
end
