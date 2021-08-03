export default class StringFinder {
    private newlinesRegex;
    constructor(needles: string[] | RegExp);
    private convertToPipedExpression(needles);
    findAll(haystack: string): {
        index: number;
        text: string;
    }[];
}
