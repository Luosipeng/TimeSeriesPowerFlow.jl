# PowerFlow 可视化示例

本文档提供了 PowerFlow 库中可视化功能的详细使用示例，展示了如何创建各种类型的图表和可视化效果，以便更好地理解和分析电力系统的计算结果。

## 基础可视化

### 系统拓扑图

创建电力系统的拓扑图：

```julia
using PowerFlow
using TimeSeriesPowerFlow
using Utils
using Plots
using GraphPlot
using LightGraphs

# 加载标准测试系统
case = load_case("case14.m")

# 创建系统拓扑图
plot_system_topology(case)

# 保存图像
savefig("system_topology.png")
```

### 自定义拓扑图

创建带有自定义布局和样式的拓扑图：

```julia
# 创建自定义拓扑图
plot_system_topology(case, 
    node_size = [degree(case, i) * 3 for i in 1:size(case.bus, 1)],  # 基于节点度的大小
    node_color = [:red for i in 1:size(case.bus, 1)],                # 自定义节点颜色
    edge_width = [case.branch[i, BR_R] * 10 + 1 for i in 1:size(case.branch, 1)],  # 基于电阻的边宽度
    layout = spring_layout,                                          # 使用弹簧布局算法
    labels = ["B$(Int(case.bus[i, BUS_I]))" for i in 1:size(case.bus, 1)]  # 添加标签
)

# 保存图像
savefig("custom_topology.png")
```

### 带有电气参数的拓扑图

创建显示电气参数的拓扑图：

```julia
# 执行潮流计算
results = runpf(case)

# 创建带有电气参数的拓扑图
plot_system_topology_with_results(case, results,
    show_voltage = true,     # 显示电压值
    show_power_flow = true,  # 显示功率流
    show_losses = true,      # 显示损耗
    edge_color_by_loading = true  # 根据负载率着色
)

# 保存图像
savefig("topology_with_results.png")
```

## 电压分析可视化

### 电压分布图

创建系统电压分布图：

```julia
# 执行潮流计算
results = runpf(case)

# 提取电压值
vm = results.bus[:, VM]
bus_ids = results.bus[:, BUS_I]

# 创建电压分布条形图
bar(bus_ids, vm,
    title = "母线电压分布",
    xlabel = "母线编号",
    ylabel = "电压幅值(pu)",
    legend = false,
    ylims = (0.9, 1.1),
    color = :blue,
    linewidth = 2,
    grid = true
)

# 添加参考线
hline!([1.0], color = :red, linestyle = :dash, label = "标称电压")
hline!([0.95, 1.05], color = :orange, linestyle = :dot, label = "限制")

# 保存图像
savefig("voltage_distribution.png")
```

### 电压等高线图

创建系统电压等高线图：

```julia
# 执行潮流计算
results = runpf(case)

# 获取母线位置信息
# 注意：实际应用中需要有母线的地理坐标
# 这里我们使用一个示例坐标
bus_x = [0, 1, 2, 3, 2, 1, 0, 1, 2, 3, 2, 1, 2, 3]
bus_y = [0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4]

# 提取电压值
vm = results.bus[:, VM]

# 创建网格
n = 100
x_grid = range(minimum(bus_x) - 0.5, maximum(bus_x) + 0.5, length = n)
y_grid = range(minimum(bus_y) - 0.5, maximum(bus_y) + 0.5, length = n)
z_grid = zeros(n, n)

# 使用反距离加权插值
for i in 1:n
    for j in 1:n
        weights = 1.0 ./ sqrt.((x_grid[i] .- bus_x).^2 + (y_grid[j] .- bus_y).^2)
        z_grid[j, i] = sum(weights .* vm) / sum(weights)
    end
end

# 创建等高线图
contour(x_grid, y_grid, z_grid,
    title = "电压等高线图",
    xlabel = "X坐标",
    ylabel = "Y坐标",
    fill = true,
    levels = 15,
    c = :viridis
)

# 添加母线位置
scatter!(bus_x, bus_y,
    marker = :circle,
    markersize = 5,
    color = :black,
    label = "母线位置"
)

# 添加母线标签
for i in 1:length(bus_x)
    annotate!(bus_x[i], bus_y[i], text("B$(Int(case.bus[i, BUS_I]))", :black, 8))
end

# 保存图像
savefig("voltage_contour.png")
```

### 3D电压曲面图

创建3D电压曲面图：

```julia
# 创建3D电压曲面图
surface(x_grid, y_grid, z_grid,
    title = "3D电压曲面图",
    xlabel = "X坐标",
    ylabel = "Y坐标",
    zlabel = "电压幅值(pu)",
    color = :viridis,
    camera = (30, 45)
)

# 添加母线位置
scatter!(bus_x, bus_y, vm,
    marker = :circle,
    markersize = 3,
    color = :red,
    label = "母线"
)

# 保存图像
savefig("voltage_surface_3d.png")
```

## 功率流可视化

### 功率流分布图

创建功率流分布图：

```julia
# 提取支路功率流
branch_flow = abs.(results.branch[:, PF])
branch_from = results.branch[:, F_BUS]
branch_to = results.branch[:, T_BUS]

# 创建功率流分布条形图
bar(1:length(branch_flow), branch_flow,
    title = "支路功率流分布",
    xlabel = "支路编号",
    ylabel = "功率流大小(MW)",
    legend = false,
    color = :green,
    linewidth = 2,
    grid = true
)

# 添加支路标签
xticks!(1:length(branch_flow), ["$(Int(branch_from[i]))-$(Int(branch_to[i]))" for i in 1:length(branch_flow)])
xtickfont(8, rotation = 45)

# 保存图像
savefig("power_flow_distribution.png")
```

### 支路负载率图

创建支路负载率图：

```julia
# 计算支路负载率
branch_rate = results.branch[:, RATE_A]
loading_rate = zeros(length(branch_flow))

for i in 1:length(branch_flow)
    if branch_rate[i] > 0
        loading_rate[i] = branch_flow[i] / branch_rate[i] * 100
    else
        loading_rate[i] = 0  # 对于没有额定值的支路
    end
end

# 创建负载率条形图
bar(1:length(loading_rate), loading_rate,
    title = "支路负载率",
    xlabel = "支路编号",
    ylabel = "负载率(%)",
    legend = false,
    color = [loading_rate[i] > 80 ? :red : (loading_rate[i] > 50 ? :orange : :green) for i in 1:length(loading_rate)],
    linewidth = 2,
    grid = true
)

# 添加支路标签
xticks!(1:length(loading_rate), ["$(Int(branch_from[i]))-$(Int(branch_to[i]))" for i in 1:length(loading_rate)])
xtickfont(8, rotation = 45)

# 添加警告线
hline!([80], color = :red, linestyle = :dash, label = "警告阈值(80%)")
hline!([100], color = :black, linestyle = :dash, label = "额定容量(100%)")

# 保存图像
savefig("branch_loading_rate.png")
```

### 功率流箭头图

创建带有功率流方向的箭头图：

```julia
# 创建带有功率流方向的箭头图
plot_power_flow_arrows(case, results,
    arrow_scale = 0.2,           # 箭头大小比例
    flow_color_by_value = true,  # 根据功率流大小着色
    show_values = true           # 显示功率流数值
)

# 保存图像
savefig("power_flow_arrows.png")
```

## 时间序列可视化

### 基本时间序列图

创建基本的时间序列图：

```julia
# 执行时间序列潮流计算
load_file = "load_data.xlsx"
ts_results = runtdpf(case, load_file)

# 提取时间序列数据
time_points = ts_results.time_points
total_load = zeros(length(time_points))
total_gen = zeros(length(time_points))
total_loss = zeros(length(time_points))
min_voltage = zeros(length(time_points))

for t in 1:length(time_points)
    if ts_results.success[t]
        result = ts_results.bus_results[t]
        total_load[t] = sum(result.bus[:, PD])
        total_gen[t] = sum(result.gen[:, PG])
        total_loss[t] = sum(result.branch[:, PL])
        min_voltage[t] = minimum(result.bus[:, VM])
    end
end

# 创建时间序列图
plot(time_points, [total_load total_gen total_loss],
    title = "系统时间序列",
    xlabel = "时间(小时)",
    ylabel = "功率(MW)",
    label = ["总负荷" "总发电" "总损耗"],
    marker = [:circle :square :diamond],
    linewidth = 2,
    grid = true,
    legend = :topleft
)

# 保存图像
savefig("time_series_basic.png")
```

### 多轴时间序列图

创建带有多个Y轴的时间序列图：

```julia
# 创建带有多个Y轴的时间序列图
p1 = plot(time_points, total_load,
    title = "系统负荷和电压时间序列",
    xlabel = "时间(小时)",
    ylabel = "负荷(MW)",
    color = :blue,
    marker = :circle,
    linewidth = 2,
    label = "总负荷",
    legend = :topleft
)

# 添加第二个Y轴
plot!(twinx(), time_points, min_voltage,
    ylabel = "电压(pu)",
    color = :red,
    marker = :square,
    linewidth = 2,
    label = "最低电压",
    legend = :topright
)

# 保存图像
savefig("time_series_dual_axis.png")
```

### 热力图

创建时间序列热力图：

```julia
# 提取所有母线的电压时间序列
voltage_ts = zeros(length(time_points), size(case.bus, 1))

for t in 1:length(time_points)
    if ts_results.success[t]
        voltage_ts[t, :] = ts_results.bus_results[t].bus[:, VM]
    end
end

# 创建电压热力图
heatmap(1:size(case.bus, 1), time_points, voltage_ts,
    title = "母线电压热力图",
    xlabel = "母线编号",
    ylabel = "时间(小时)",
    color = :viridis,
    colorbar_title = "电压(pu)"
)

# 保存图像
savefig("voltage_heatmap.png")

# 提取所有支路的负载率时间序列
loading_ts = zeros(length(time_points), size(case.branch, 1))

for t in 1:length(time_points)
    if ts_results.success[t]
        for i in 1:size(case.branch, 1)
            flow = abs(ts_results.bus_results[t].branch[i, PF])
            rate = case.branch[i, RATE_A]
            if rate > 0
                loading_ts[t, i] = flow / rate * 100
            end
        end
    end
end

# 创建支路负载率热力图
heatmap(1:size(case.branch, 1), time_points, loading_ts,
    title = "支路负载率热力图",
    xlabel = "支路编号",
    ylabel = "时间(小时)",
    color = :thermal,
    colorbar_title = "负载率(%)"
)

# 保存图像
savefig("loading_heatmap.png")
```

## 高级可视化

### 桑基图

创建电力流向的桑基图：

```julia
using SankeyPlots

# 执行潮流计算
results = runpf(case)

# 准备桑基图数据
# 提取发电机、负荷和支路数据
gen_buses = results.gen[:, GEN_BUS]
gen_output = results.gen[:, PG]

# 创建节点列表
nodes = []
# 添加发电节点
for i in 1:length(gen_buses)
    push!(nodes, "G$(i)_$(Int(gen_buses[i]))")
end
# 添加母线节点
for i in 1:size(results.bus, 1)
    push!(nodes, "B$(Int(results.bus[i, BUS_I]))")
end
# 添加负荷节点
load_buses = findall(results.bus[:, PD] .> 0)
for i in load_buses
    push!(nodes, "L$(Int(results.bus[i, BUS_I]))")
end

# 创建链接列表
links = []
link_values = []

# 发电机到母线的链接
for i in 1:length(gen_buses)
    gen_idx = "G$(i)_$(Int(gen_buses[i]))"
    bus_idx = "B$(Int(gen_buses[i]))"
    push!(links, (gen_idx, bus_idx))
    push!(link_values, gen_output[i])
end

# 母线之间的链接(支路功率流)
for i in 1:size(results.branch, 1)
    from_bus = Int(results.branch[i, F_BUS])
    to_bus = Int(results.branch[i, T_BUS])
    power_flow = results.branch[i, PF]
    
    if power_flow > 0
        push!(links, ("B$(from_bus)", "B$(to_bus)"))
        push!(link_values, abs(power_flow))
    else
        push!(links, ("B$(to_bus)", "B$(from_bus)"))
        push!(link_values, abs(power_flow))
    end
end

# 母线到负荷的链接
for i in load_buses
    bus_id = Int(results.bus[i, BUS_I])
    load = results.bus[i, PD]
    push!(links, ("B$(bus_id)", "L$(bus_id)"))
    push!(link_values, load)
end

# 创建桑基图
sankey(links, link_values,
    title = "电力流向桑基图",
    node_labels = nodes,
    node_colors = [i <= length(gen_buses) ? :green : 
                  (i <= length(gen_buses) + size(results.bus, 1) ? :blue : :red) 
                  for i in 1:length(nodes)],
    link_color = :gray
)

# 保存图像
savefig("power_flow_sankey.png")
```

### 雷达图

创建系统性能雷达图：

```julia
# 执行潮流计算
results = runpf(case)

# 计算系统性能指标
voltage_deviation = std(results.bus[:, VM])
voltage_range = maximum(results.bus[:, VM]) - minimum(results.bus[:, VM])
total_loss = sum(results.branch[:, PL])
loss_percent = total_loss / sum(results.gen[:, PG]) * 100
max_loading = maximum(abs.(results.branch[:, PF]) ./ results.branch[:, RATE_A] .* 100)
gen_reserve = sum(results.gen[:, PMAX] - results.gen[:, PG]) / sum(results.gen[:, PMAX]) * 100

# 归一化指标(越低越好)
max_values = [0.05, 0.2, 50, 10, 100, 100]  # 各指标的最大参考值
indicators = [voltage_deviation, voltage_range, total_loss, loss_percent, max_loading, 100 - gen_reserve]
normalized = min.(indicators ./ max_values, 1)

# 创建雷达图
labels = ["电压偏差", "电压范围", "总损耗(MW)", "损耗率(%)", "最大负载率(%)", "发电预留(%)"]
radar_plot(normalized, labels,
    title = "系统性能雷达图",
    fill = true,
    color = :blue,
    linewidth = 2,
    grid = true
)

# 保存图像
savefig("system_performance_radar.png")
```

### 地理信息系统(GIS)可视化

使用地理信息创建系统可视化：

```julia
using GeoPlots

# 假设我们有母线的地理坐标
# 实际应用中应从GIS数据库获取
bus_lat = [40.0, 40.1, 40.2, 40.1, 39.9, 39.8, 39.7, 39.8, 39.9, 40.0, 40.1, 40.2, 40.3, 40.4]
bus_lon = [-74.0, -74.1, -74.2, -74.3, -74.2, -74.1, -74.0, -73.9, -73.8, -73.7, -73.8, -73.9, -74.0, -74.1]

# 创建地理地图
geoplot(bus_lat, bus_lon,
    title = "电力系统地理分布",
    marker = :circle,
    markersize = 8,
    color = :blue,
    label = "母线位置"
)

# 添加支路连接
for i in 1:size(case.branch, 1)
    from_idx = findfirst(x -> x == case.branch[i, F_BUS], case.bus[:, BUS_I])
    to_idx = findfirst(x -> x == case.branch[i, T_BUS], case.bus[:, BUS_I])
    
    if from_idx !== nothing && to_idx !== nothing
        plot!([bus_lat[from_idx], bus_lat[to_idx]], [bus_lon[from_idx], bus_lon[to_idx]],
            color = :gray,
            linewidth = 2,
            label = nothing
        )
    end
end

# 添加母线标签
for i in 1:length(bus_lat)
    annotate!(bus_lat[i], bus_lon[i], text("B$(Int(case.bus[i, BUS_I]))", :black, 8))
end

# 保存图像
savefig("geographic_system_map.png")
```

## 动态可视化

### 时间序列动画

创建电压变化的动画：

```julia
using Plots

# 执行时间序列潮流计算
ts_results = runtdpf(case, load_file)

# 创建电压分布动画
anim = @animate for t in 1:length(ts_results.time_points)
    if ts_results.success[t]
        result = ts_results.bus_results[t]
        bar(result.bus[:, BUS_I], result.bus[:, VM],
            title = "电压分布 - 时间点 $(ts_results.time_points[t])",
            xlabel = "母线编号",
            ylabel = "电压幅值(pu)",
            ylims = (0.9, 1.1),
            legend = false,
            color = :blue
        )
        # 添加参考线
        hline!([1.0], color = :red, linestyle = :dash, label = "标称电压")
        hline!([0.95, 1.05], color = :orange, linestyle = :dot, label = "限制")
    end
end

# 保存动画
gif(anim, "voltage_animation.gif", fps = 2)
```

### 功率流动画

创建功率流变化的动画：

```julia
# 创建功率流动画
anim = @animate for t in 1:length(ts_results.time_points)
    if ts_results.success[t]
        result = ts_results.bus_results[t]
        branch_flow = abs.(result.branch[:, PF])
        
        bar(1:length(branch_flow), branch_flow,
            title = "支路功率流 - 时间点 $(ts_results.time_points[t])",
            xlabel = "支路编号",
            ylabel = "功率流大小(MW)",
            legend = false,
            color = :green
        )
    end
end

# 保存动画
gif(anim, "power_flow_animation.gif", fps = 2)
```

### 交互式拓扑图

创建交互式系统拓扑图：

```julia
using PlotlyJS

# 准备节点数据
node_x = bus_x
node_y = bus_y
node_text = ["母线 $(Int(case.bus[i, BUS_I]))<br>电压: $(round(results.bus[i, VM], digits=3)) pu<br>负荷: $(round(case.bus[i, PD], digits=1)) MW" for i in 1:length(bus_x)]
node_size = [case.bus[i, PD] > 0 ? 10 + case.bus[i, PD] / 2 : 10 for i in 1:length(bus_x)]
node_color = [case.bus[i, BUS_TYPE] == 3 ? "red" : (case.bus[i, BUS_TYPE] == 2 ? "blue" : "gray") for i in 1:length(bus_x)]

# 准备边数据
edge_x = []
edge_y = []
edge_text = []

for i in 1:size(case.branch, 1)
    from_idx = findfirst(x -> x == case.branch[i, F_BUS], case.bus[:, BUS_I])
    to_idx = findfirst(x -> x == case.branch[i, T_BUS], case.bus[:, BUS_I])
    
    if from_idx !== nothing && to_idx !== nothing
        push!(edge_x, bus_x[from_idx], bus_x[to_idx], nothing)
        push!(edge_y, bus_y[from_idx], bus_y[to_idx], nothing)
        
        flow = abs(results.branch[i, PF])
        loss = results.branch[i, PL]
        text = "支路 $(Int(case.branch[i, F_BUS]))-$(Int(case.branch[i, T_BUS]))<br>功率流: $(round(flow, digits=2)) MW<br>损耗: $(round(loss, digits=2)) MW"
        push!(edge_text, text, text, nothing)
    end
end

# 创建边迹线
edge_trace = scatter(
    x = edge_x, y = edge_y,
    line = attr(width = 1.5, color = "#888"),
    hoverinfo = "text",
    text = edge_text,
    mode = "lines"
)

# 创建节点迹线
node_trace = scatter(
    x = node_x, y = node_y,
    mode = "markers",
    hoverinfo = "text",
    text = node_text,
    marker = attr(
        showscale = true,
        colorscale = "YlGnBu",
        color = node_color,
        size = node_size,
        line = attr(width = 2, color = "black")
    )
)

# 创建布局
layout = Layout(
    title = "交互式电力系统拓扑图",
    showlegend = false,
    hovermode = "closest",
    margin = attr(b = 20, l = 5, r = 5, t = 40),
    xaxis = attr(showgrid = false, zeroline = false, showticklabels = false),
    yaxis = attr(showgrid = false, zeroline = false, showticklabels = false)
)

# 创建图表
p = plot([edge_trace, node_trace], layout)

# 保存为HTML文件
savefig(p, "interactive_topology.html")
```

## 违例分析可视化

### 电压违例可视化

可视化电压违例：

```julia
# 执行时间序列潮流计算
ts_results = runtdpf(case, load_file)

# 设置电压限制
v_min = 0.95
v_max = 1.05

# 提取电压数据
voltage_data = zeros(length(ts_results.time_points), size(case.bus, 1))
for t in 1:length(ts_results.time_points)
    if ts_results.success[t]
        voltage_data[t, :] = ts_results.bus_results[t].bus[:, VM]
    end
end

# 创建电压违例热力图
voltage_violations = (voltage_data .< v_min) .| (voltage_data .> v_max)
violation_mask = zeros(size(voltage_data))
violation_mask[voltage_violations] .= 1

# 使用掩码创建热力图
heatmap(1:size(case.bus, 1), ts_results.time_points, violation_mask,
    title = "电压违例热力图",
    xlabel = "母线编号",
    ylabel = "时间(小时)",
    color = cgrad([:white, :red]),
    colorbar_title = "违例状态"
)

# 保存图像
savefig("voltage_violations_heatmap.png")

# 创建违例计数条形图
violation_count = sum(voltage_violations, dims=1)[:]
bar(1:size(case.bus, 1), violation_count,
    title = "电压违例计数",
    xlabel = "母线编号",
    ylabel = "违例时间点数量",
    color = :red,
    legend = false
)

# 保存图像
savefig("voltage_violations_count.png")
```

### 支路过载可视化

可视化支路过载情况：

```julia
# 设置负载率阈值
loading_threshold = 80  # 80%

# 提取支路负载率数据
loading_data = zeros(length(ts_results.time_points), size(case.branch, 1))
for t in 1:length(ts_results.time_points)
    if ts_results.success[t]
        for i in 1:size(case.branch, 1)
            flow = abs(ts_results.bus_results[t].branch[i, PF])
            rate = case.branch[i, RATE_A]
            if rate > 0
                loading_data[t, i] = flow / rate * 100
            end
        end
    end
end

# 创建支路过载热力图
loading_violations = loading_data .> loading_threshold
violation_mask = zeros(size(loading_data))
violation_mask[loading_violations] .= 1

# 使用掩码创建热力图
heatmap(1:size(case.branch, 1), ts_results.time_points, violation_mask,
    title = "支路过载热力图(>$(loading_threshold)%)",
    xlabel = "支路编号",
    ylabel = "时间(小时)",
    color = cgrad([:white, :orange]),
    colorbar_title = "违例状态"
)

# 保存图像
savefig("loading_violations_heatmap.png")

# 创建违例计数条形图
violation_count = sum(loading_violations, dims=1)[:]
bar(1:size(case.branch, 1), violation_count,
    title = "支路过载计数",
    xlabel = "支路编号",
    ylabel = "违例时间点数量",
    color = :orange,
    legend = false
)

