function shuffledList=fisheryates(inputList)
% Fisher Yates (aka Knuth) shuffle
% Make a copy of the input list to avoid modifying the original list
shuffledList = inputList;
% Start from the end of the list and swap each element with a randomly chosen element
n = length(shuffledList);
for i = n:-1:2
    j = randi(i-1);  % Generate a random index between 1 and i (inclusive)
    % Swap elements
    temp = shuffledList(i);
    shuffledList(i) = shuffledList(j);
    shuffledList(j) = temp;
end
end