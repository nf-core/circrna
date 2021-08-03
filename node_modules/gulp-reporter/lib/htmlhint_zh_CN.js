'use strict';

const messages = new Map([
	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/alt-require.js
	['An alt attribute must be present on <img> elements.', '<img>元素必须要有alt属性。'],
	[/^The alt attribute of (.+?) must have a value\.$/, '$1必须要有alt属性。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-lowercase.js
	[/^The attribute name of (.+?) must be in lowercase\.$/, '属性名$1必须小写。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-no-duplication.js
	[/^Duplicate of attribute name (.+?) was found\.$/, '发现重复的属性$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-unsafe-chars.js
	[/^The value of attribute (.+?) cannot contain an unsafe char (.+?)\.$/, '属性$1的值中不能包含不安全字符$2。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-value-double-quotes.js
	[/^The value of attribute (.+?) must be in double quotes\.$/, '属性$1的值必须放在双引号中。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/attr-value-not-empty.js
	[/^The attribute (.+?) must have a value\.$/, '属性$1必须有值。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/doctype-first.js
	['Doctype must be declared first.', 'doctype必须首先声明。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/doctype-html5.js
	['Invalid doctype. Use: "<!DOCTYPE html>"', '无效的doctype，请使用：“<!DOCTYPE html>”。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/href-abs-or-rel.js
	[/^The value of the href attribute (.+?) must be absolute\.$/, '属性href的值$1必须是绝对路径。'],
	[/^The value of the href attribute (.+?) must be relative\.$/, '属性href的值$1必须是相对路径。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/id-class-ad-disabled.js
	[/^The value of attribute (.+?) cannot use the ad keyword\.$/, '属性$1的值不能包含关键字“ad”。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/id-class-value.js
	['The id and class attribute values must be in lowercase and split by an underscore.', 'id和class属性的值必须小写，并由下划线分割。'],
	['The id and class attribute values must be in lowercase and split by a dash.', 'id和class属性的值必须小写，并由连字符分割。'],
	['The id and class attribute values must meet the camelCase style.', 'id和class属性的值必须符合camelCase风格。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/id-unique.js
	[/^The id value (.+?) must be unique\.$/, 'id的值$1必须是唯一的。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/inline-script-disabled.js
	[/^Inline script (.+?) cannot be used\.$/, '不能使用行内脚本$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/inline-style-disabled.js
	[/^Inline style (.+?) cannot be used\.$/, '不能使用行内样式$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/space-tab-mixed-disabled.js
	[/^Please use space for indentation and keep (.+?) length\.$/, '请使用$1个空格作为缩进。'],
	['Please use space for indentation.', '请使用空格作为缩进。'],
	['Please use tab for indentation.', '请使用tab作为缩进。'],
	['Do not mix tabs and spaces for indentation.', '不要混用tab和空格作为缩进。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/spec-char-escape.js
	[/^Special characters must be escaped ?: (.+?)\.$/, '必须转义特殊字符：$1。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/src-not-empty.js
	[/^The attribute (.+?) of the tag (.+?) must have a value\.$/, '$1标签的属性$2必须有值。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/style-disabled.js
	['The <style> tag cannot be used.', '不能使用<style>标签。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/tag-pair.js
	[/^Tag must be paired, missing: (.+?), (?:start|open) tag match failed (.+?) on line (.+?)\.$/, '标签必须匹配，缺失：$1，在第$3行匹配开始标签$2失败'],
	[/^Tag must be paired, no start tag: (.+?)$/, '标签必须匹配，没有开始标签：$1'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/tag-self-close.js
	[/^The empty tag ?: (.+?) must be self closed\.$/, '空标签：$1必须自关闭'],

	[/^The html element name of (.+?) must be in lowercase\.$/, 'html元素$1的名字必须小写。'],

	// https://github.com/yaniswang/HTMLHint/blob/master/src/rules/title-require.js
	['<title></title> must not be empty.', '<title></title> 必须不为空。'],
	['<title> must be present in <head> tag.', '<title>必须存在于<head>标签。'],
]);
module.exports = function (error) {
	for (const [en, han] of messages) {
		if (en instanceof RegExp) {
			if (en.test(error.message)) {
				error.message = error.message.replace(en, han);
				break;
			}
		} else if (error.message === en) {
			error.message = han;
			break;
		}
	}
	return error;
};
