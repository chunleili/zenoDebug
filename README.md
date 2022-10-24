# 用于zeno debug的模板

zeno由于采用tcp的设计导致调试器都无法使用，这会给debug带来极大痛苦。

因此我们提炼出implementation代码，让其与zeno彻底分离。从而使用调试器。

使用方法是：
1. 在zeno中，建立节点。该节点只是传入数据的，比如pos。并且传的就是单纯的std::vector。然后传出数据，用于可视化的。
2. 把所有具体实现的代码抽离出来，放到类Impl中。相应的代码拷贝到Impl.cpp和Impl.h。这里不要用任何除了vec以外的zeno相关内容。
3. 在zenoDebug（本项目）中运行调试器。并进行调试。
4. 将Impl类中的调试完毕的代码拷贝回zeno中。