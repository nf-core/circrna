'use strict';
/**
 * 将错误等级的字符串转化为排序所用的权重
 * @param  {String} severity 错误等级
 * @return {Int}          排序权重
 */
function severity2num (severity) {
	if (!severity) {
		severity = 0;
	} else if (severity === 'error') {
		severity = 1;
	} else if (severity === 'warn') {
		severity = 2;
	} else if (severity === 'info') {
		severity = 3;
	} else {
		severity = 4;
	}
	return severity;
}

/**
 * 将错误按照severity>lineNumber>columnNumber这三个维度的优先级排序时，所用的的对象比较函数。
 *
 * @param {Error} err1 要比较的对象
 * @param {Error} err2 要比较的对象
 * @returns {int} 比较结果，正数则err1大，负数err2大，0相等
 */
function compareDefault (err1, err2) {
	if (err1.fileName !== err2.fileName) {
		return (err1.fileName || '').localeCompare(err2.fileName || '');
	} else if (err1.demote !== err2.demote) {
		return err1.demote ? 1 : -1;
	} else if (err1.severity !== err2.severity) {
		return severity2num(err1.severity) - severity2num(err2.severity);
	} else if (err1.lineNumber === err2.lineNumber) {
		return (err1.columnNumber || 0) - (err2.columnNumber || 0);
	} else {
		return (err1.lineNumber || 0) - (err2.lineNumber || 0);
	}
}

/**
 * 对错误信息排序的函数
 * @param  {Array} errors 待排序的错误数组
 * @param  {Function} compareFunction 可选。用来指定按某种顺序进行排列的函数。如果省略，元素按照severity>lineNumber>columnNumber这三个维度的优先级排序。
 * @return {Array}        排好序的错误数组
 */
function sortError (errors, compareFunction) {
	compareFunction = compareFunction || compareDefault;
	return errors.sort(compareFunction);
}

module.exports = sortError;
