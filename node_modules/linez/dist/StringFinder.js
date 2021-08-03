"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var StringFinder = (function () {
    function StringFinder(needles) {
        if (needles instanceof RegExp) {
            if (!needles.global) {
                throw new Error('StringFinder regular expression must have a global flag');
            }
            this.newlinesRegex = needles;
            return;
        }
        if (needles instanceof Array) {
            this.newlinesRegex = this.convertToPipedExpression(needles);
            return;
        }
        throw new Error('Unexpected type in StringFinder constructor argument: ' + typeof needles);
    }
    StringFinder.prototype.convertToPipedExpression = function (needles) {
        needles = needles.map(function (needle) {
            return needle.replace('\\', '\\\\');
        });
        return new RegExp('(' + needles.join('|') + ')', 'g');
    };
    StringFinder.prototype.findAll = function (haystack) {
        var matches = [];
        while (true) {
            var match = this.newlinesRegex.exec(haystack);
            if (!match) {
                break;
            }
            matches.push({
                index: match.index,
                text: match[0]
            });
        }
        return matches;
    };
    return StringFinder;
}());
exports.default = StringFinder;
//# sourceMappingURL=StringFinder.js.map