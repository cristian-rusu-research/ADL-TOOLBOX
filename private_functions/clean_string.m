function y = clean_string(y)
y(y == ' ') = [];
y = upper(y);
