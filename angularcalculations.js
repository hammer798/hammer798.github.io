// JavaScript source code
let app = angular.module('phyApp', []);

app.controller('myCtrl', function ($scope) {
    $scope.functions = [
        JSON.parse('{"name": "kin0","value":"", "orderOps":"=+*", "vars":["v", "v0", "a", "t"]}'),
        JSON.parse('{ "name": "notpotato", "value": "5"}')
    ];

    $scope.getValues = function (f) {
        if (f == null || f.vars == null || f.vars.length < 1)
            return;
        (f.vars).forEach(function (v) {
            v.value = getData(v);
        })
    }

});

//pulls data from relevant html input element
function getData(id) {
    id = id.toString();

    let docElement = document.getElementById(id);

    let data = null;
    if(docElement != null && docElement.value != null)
        data = parseFloat(docElement.value);
    return data;
}

