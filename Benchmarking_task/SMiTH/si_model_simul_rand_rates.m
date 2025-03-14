% This script simulates random transmission trees under the SI model with
% variable transmission rates. Specific details are described in the paper

clear;
nTests = 500;

for t = 1:nTests
    foldIn = ['Favites_new' filesep 'FAVITES_output_expSI_contemp_T200_N100_E1_' int2str(t)];
    foldOut = ['Favites_new_revision' filesep 'TransNet_' int2str(t)];
    cn_file_in = [foldIn filesep 'contact_network.txt'];
    cn_file_out = [foldOut filesep 'contact_network.txt'];
    tn_file_out = [foldOut filesep 'trans_network_randrate.txt'];
    config_file_in = [foldIn filesep 'CONFIG.json'];
    config_file_out = [foldOut filesep 'CONFIG_fixedNet.json'];
    log_file = [foldOut filesep 'log.csv'];
    
    a = 0.002;
    b = 0.05;
    maxTime = 100;

    if ~isfolder(foldOut)
        mkdir(foldOut);
    end

    toTry = true;
    while toTry
        [trans_net,trans_rates] = si_model_simul(cn_file_in, tn_file_out, a, b,maxTime);
        if ismember(length(trans_net),5:30) 
            toTry = false;
        end
    end


    fid = fopen(config_file_in, 'rt', 'n', 'UTF-8'); 
    str = fscanf(fid, '%c'); 
    fclose(fid);
    str = regexprep(str, '''', '"');
    str = regexprep(str, '^\xEF\xBB\xBF', '');
    str = strtrim(str);


    m = length(str);
    str(end-18:end-1) = '';
    data = jsondecode(str);
    data.TransmissionNodeSample = "TransmissionFile";
    data.transmission_network_file = tn_file_out;
    data.EndCriteria = "TransmissionFile";
    data.SeedSelection = "TransmissionFile";
    data.TransmissionTimeSample = "TransmissionFile";
    data.ContactNetworkGenerator = "File";
    data.contact_network_file = cn_file_out;
    data.NodeEvolution = 'BirthDeath';
    data.bd_birth = 0.1;
    data.bd_birth_sd = 0;
    data.bd_death = 0;
    data.bd_death_sd = 0;

    c= data.seed_seqs;
    c = c{1};
    ind = true(1,length(c));
    for i = 1:length(c)
        if (c(i)~='A') && (c(i)~='C') && (c(i)~='T') && (c(i)~='G')
            ind(i) = 0;
        end
    end
    c = c(ind);
    data.seed_seqs={c};

    outData = jsonencode(data);
    fid = fopen(config_file_out, 'w', 'n', 'UTF-8');
    fwrite(fid, outData, 'char');
    fclose(fid);
    copyfile(cn_file_in, cn_file_out);
    writematrix(trans_rates,log_file);
end


function [transmission_network,transmission_rates] = si_model_simul(contact_network_file, output_file, a, b,maxTime)
    % Load the contact network from the file
    [nodes, edges] = load_contact_network(contact_network_file);

    % Generate random transmission rates for edges
    num_edges = size(edges, 1);
    transmission_rates = a + (b - a) * rand(num_edges, 1);

    % Compute degrees for each node
    degrees = compute_degrees(nodes, edges);

    % Select the initial infected node based on degree-proportional probabilities
    initial_node = select_initial_node(degrees);

    % Simulate the SI model
    transmission_network = simulate_si(nodes, edges, transmission_rates, initial_node, maxTime);
    

    % Save the transmission network to the output file
    save_transmission_network(output_file, transmission_network,initial_node);
end

function [nodes, edges] = load_contact_network(filename)
    fileID = fopen(filename, 'r');
    lines = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    
    lines = lines{1};
    nodes = [];
    edges = [];
    
    for i = 1:length(lines)
        line = strsplit(lines{i}, '\t');
        if strcmp(line{1}, 'NODE')
            nodes(end + 1) = str2double(line{2});
        elseif strcmp(line{1}, 'EDGE')
            edges(end + 1, :) = [str2double(line{2}), str2double(line{3})];
        end
    end
end

function degrees = compute_degrees(nodes, edges)
    degrees = zeros(max(nodes) + 1, 1);
    for i = 1:size(edges, 1)
        degrees(edges(i, 1) + 1) = degrees(edges(i, 1) + 1) + 1;
        degrees(edges(i, 2) + 1) = degrees(edges(i, 2) + 1) + 1;
    end
end

function initial_node = select_initial_node(degrees)
    probabilities = degrees / sum(degrees);
    initial_node = find(rand < cumsum(probabilities), 1) - 1;
end

function transmission_network = simulate_si(nodes, edges, rates, initial_node, maxTime)
    infected = false(max(nodes) + 1, 1);
    infected(initial_node + 1) = true;

    transmission_network = {};
    transmission_network{end + 1} = {'None', initial_node, 0};
    event_queue = struct('time', [], 'source', [], 'target', []);

    % Initialize the event queue with potential transmissions from the initial node
    for i = 1:size(edges, 1)
        if edges(i, 1) == initial_node
            time = exprnd(1 / rates(i));
            event_queue = add_event(event_queue, time, edges(i, 1), edges(i, 2));
        elseif edges(i, 2) == initial_node
            time = exprnd(1 / rates(i));
            event_queue = add_event(event_queue, time, edges(i, 2), edges(i, 1));
        end
    end

    % Process the event queue
    while ~isempty(event_queue.time)
        [time, idx] = min(event_queue.time);

        % Stop the simulation if the next event exceeds maxTime
        if time > maxTime
            break;
        end

        source = event_queue.source(idx);
        target = event_queue.target(idx);

        % Remove the processed event
        event_queue = remove_event(event_queue, idx);

        % Only process the event if the source is infected
        if infected(source + 1) && ~infected(target + 1)
            infected(target + 1) = true;
            transmission_network{end + 1} = [source, target, time];

            % Add new transmission events for the newly infected node
            for i = 1:size(edges, 1)
                if edges(i, 1) == target && ~infected(edges(i, 2) + 1)
                    time = time + exprnd(1 / rates(i));
                    event_queue = add_event(event_queue, time, edges(i, 1), edges(i, 2));
                elseif edges(i, 2) == target && ~infected(edges(i, 1) + 1)
                    time = time + exprnd(1 / rates(i));
                    event_queue = add_event(event_queue, time, edges(i, 2), edges(i, 1));
                end
            end
        end
    end
end

function event_queue = add_event(event_queue, time, source, target)
    event_queue.time = [event_queue.time; time];
    event_queue.source = [event_queue.source; source];
    event_queue.target = [event_queue.target; target];
end

function event_queue = remove_event(event_queue, idx)
    event_queue.time(idx) = [];
    event_queue.source(idx) = [];
    event_queue.target(idx) = [];
end

function save_transmission_network(filename, transmission_network,initial_node)
    fileID = fopen(filename, 'w');
    fprintf(fileID, '%s\t%d\t%f\n', 'None', initial_node, 0);
    for i = 2:numel(transmission_network)
        entry = transmission_network{i};
        fprintf(fileID, '%d\t%d\t%f\n', entry(1), entry(2), entry(3));
    end
    fclose(fileID);
end